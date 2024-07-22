import os
import datetime
import threading
import subprocess
import logging
from flask import Flask, request, jsonify
from flask_cors import CORS
import joblib
import pandas as pd
from Bio import Entrez, Medline
import importlib.util  # For dynamic module loading

# Configuration
MODEL_PATH = "interaction_classifier.pkl"
VECTORIZER_PATH = "vectorizer.pkl"
LOG_FILE = "interaction_log.txt"

logging.basicConfig(filename="error.log", level=logging.ERROR)

app = Flask(__name__)
CORS(app, resources={r"/process_drugs": {"origins": "*", "methods": ["POST", "OPTIONS"]}})

# Load models and vectorizer (using dynamic import for better organization)
def load_model_and_vectorizer():
    try:
        spec = importlib.util.spec_from_file_location("drug_model", "drug_model.py")
        drug_model_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(drug_model_module)
        model, vectorizer = drug_model_module.load_model(MODEL_PATH, VECTORIZER_PATH)
        print("Loaded pre-trained model and vectorizer.")
        return model, vectorizer, drug_model_module.predict_interactions  # Get the prediction function as well
    except (FileNotFoundError, AttributeError) as e:
        logging.error(f"Error loading models or functions: {e}")
        print(f"Error loading models or functions. Check file paths and drug_model.py")
        exit(1)

model, vectorizer, predict_interactions = load_model_and_vectorizer()  # Load outside the route

def fetch_pubmed_articles(term, retmax=100):
    """Fetch PubMed articles with improved error handling."""
    from Bio import Entrez, Medline
    Entrez.email = "sannody0@outlook.sa"  # Your email
    try:
        handle = Entrez.esearch(db="pubmed", term=term, retmax=retmax)
        record = Entrez.read(handle)
        id_list = record["IdList"]
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        records = Medline.parse(handle)
        return list(records)
    except Exception as e:
        logging.error(f"Error fetching PubMed articles: {e}")  # Log the error
        return []  # Return an empty list on error

def preprocess_articles(articles):
    """Preprocess the fetched articles."""
    preprocessed_articles = []
    for record in articles:
        title = record.get('TI', 'No title')
        abstract = record.get('AB', 'No abstract')
        if abstract != 'No abstract':  # Ensure only valid abstracts are used
            preprocessed_articles.append({'title': title, 'abstract': abstract})
    return preprocessed_articles


def predict_interactions(model, vectorizer, articles):
    df = pd.DataFrame(articles)
    X = vectorizer.transform(df['abstract'])

    # Print feature consistency check
    print(f"Number of features in vectorizer: {X.shape[1]}")
    print(f"Expected number of features: {vectorizer.get_feature_names_out().shape[0]}") 

    try:
        # Try to use predict_proba if available
        y_prob = model.predict_proba(X)[:, 1]
    except AttributeError as e:
        logging.warning(f"Model does not support predict_proba: {e}")
        y_prob = [None] * X.shape[0]  # Assign None or 0.5 as default probability

    predictions = model.predict(X)
    return predictions, y_prob

@app.route('/process_drugs', methods=['POST'])
def process_drugs():
    print("Request received")
    try:
        data = request.get_json()
        print("Data received:", data)
        
        if not data:
            raise ValueError("No data received")

        drug1, drug2 = data.get('drug1'), data.get('drug2')
        if not drug1 or not drug2:
            raise ValueError("Missing drug names in the request.")

        print("Drugs:", drug1, drug2)

        drug1, drug2 = drug1.strip(), drug2.strip()
        search_term = f"({drug1} OR {drug2}) (drug interaction OR adverse reaction OR side effect)"
        print("Search term:", search_term)
        
        articles = fetch_pubmed_articles(search_term)
        print(f"Fetched {len(articles)} articles for terms: {drug1}, {drug2}")

        if not articles:
            return jsonify({"message": "No relevant articles found."}), 404

        processed_articles = preprocess_articles(articles)
        print(f"Processed articles: {len(processed_articles)}")

        predictions, probabilities = predict_interactions(model, vectorizer, processed_articles)
        print(f"Predicted interactions for {len(processed_articles)} articles")

        pubmed_results = []
        for i, article in enumerate(processed_articles):
            pubmed_results.append({
                "title": article.get("title", ""),
                "prediction": int(predictions[i]),  # Convert to int
                "probability": probabilities[i] if probabilities[i] is not None else "N/A"

            })

        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(LOG_FILE, "a") as log_file:
            log_file.write(
                f"{timestamp} - Drugs: {drug1}, {drug2} - PubMed Results: {pubmed_results}\n"
            )

        return jsonify({"pubmed_results": pubmed_results})

    except ValueError as ve:
        logging.error(f"ValueError: {ve}")
        print(f"ValueError: {ve}")
        return jsonify({"error": str(ve)}), 400
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        print(f"Unexpected error: {e}")
        return jsonify({"error": "Internal Server Error"}), 500

def start_flask_app():
    app.run(debug=True, port=os.getenv('PORT', 8000))

def start_gunicorn():
    subprocess.Popen(["gunicorn", "process_drugs:app", "-b", f"0.0.0.0:{os.getenv('PORT', '8000')}", "--reload", "--log-level", "debug"])  

if __name__ == "__main__":
    threading.Thread(target=start_gunicorn).start()
