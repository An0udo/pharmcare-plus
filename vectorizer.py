import joblib

# Load the trained vectorizer
loaded_vectorizer = joblib.load("vectorizer.pkl")

# Example new abstract to transform
new_abstract = "A new abstract about drug safety."

# Transform the new abstract using the loaded vectorizer
tfidf_representation = loaded_vectorizer.transform([new_abstract])

# The resulting 'tfidf_representation' is a sparse matrix containing the TF-IDF representation of the 'new_abstract'
