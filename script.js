async function fetchSuggestions(input, suggestionBoxId) {
    const query = input.value;
    if (query.length < 2) {
        document.getElementById(suggestionBoxId).innerHTML = '';
        return;
    }

    const response = await fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${query}/synonyms/JSON`);
    const data = await response.json();
    let suggestions = [];

    if (data.InformationList && data.InformationList.Information) {
        suggestions = data.InformationList.Information.flatMap(info => info.Synonym || []);
    }

    let scientificNames = [...new Set(suggestions)].slice(0, 5);

    let suggestionsHtml = '';
    scientificNames.forEach(suggestion => {
        if (suggestion.length <= 30) {
            suggestionsHtml += `<div onclick="selectSuggestion('${suggestion}', '${input.id}', '${suggestionBoxId}')">${suggestion}</div>`;
        }
    });
    document.getElementById(suggestionBoxId).innerHTML = suggestionsHtml;
}

function selectSuggestion(suggestion, inputId, suggestionBoxId) {
    document.getElementById(inputId).value = suggestion;
    document.getElementById(suggestionBoxId).innerHTML = '';
}

function checkInteraction() {
    const drug1 = document.getElementById('drug1').value;
    const drug2 = document.getElementById('drug2').value;
    sendDrugsToPython(drug1, drug2);
}

function clearResult() {
    document.getElementById('warningContainer').style.display = 'block';
    document.getElementById('resultContainer').style.display = 'none';
    document.getElementById('resultText').innerText = '';
}

async function sendDrugsToPython(drug1, drug2) {
    try {
        const response = await fetch('http://127.0.0.1:8009/process_drugs', { 
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({ drug1, drug2 }),
        });

        if (response.ok) {
            const result = await response.json();
            displayResults(result.pubmed_results);
        } else {
            handleError(response);
        }
    } catch (error) { 
        console.error("Error:", error);
        document.getElementById('resultText').innerHTML = '<p>Failed to connect to the server or an unexpected error occurred. Please try again.</p>';
        document.getElementById('warningContainer').style.display = 'none';
        document.getElementById('resultContainer').style.display = 'block';
    }
}

function displayResults(results) {
    let highestProbability = 0;
    let summaryText = '';

    results.forEach(result => {
        if (result.probability > highestProbability) {
            highestProbability = result.probability;
        }

        summaryText += `${result.title}: ${result.probability.toFixed(4)}. `;
    });

    // Limit the summary to 250 words
    const summaryWords = summaryText.split(' ');
    if (summaryWords.length > 250) {
        summaryText = summaryWords.slice(0, 250).join(' ') + '...';
    }

    let riskLevel = '';
    if (highestProbability >= 0.8) {
        riskLevel = 'خطر جدا';
    } else if (highestProbability >= 0.5) {
        riskLevel = 'خطر متوسط';
    } else {
        riskLevel = 'خطر ضعيف';
    }

    let resultHtml = `<p>مستوى الخطر: ${riskLevel}</p>`;
    resultHtml += `<p>ملخص الأعراض الجانبية: ${summaryText}</p>`;

    document.getElementById('resultText').innerHTML = resultHtml;
    document.getElementById('warningContainer').style.display = 'none';
    document.getElementById('resultContainer').style.display = 'block';
}

async function handleError(response) {
    if (response.status === 400) { // Bad Request (e.g., invalid input)
        const errorData = await response.json();
        document.getElementById('resultText').innerHTML = `<p>Error: ${errorData.error}</p>`;
    } else if (response.status === 404) { // Not Found (e.g., no articles found)
        document.getElementById('resultText').innerHTML = "<p>No relevant PubMed articles found.</p>";
    } else if (response.status >= 500) { // Server Error
        document.getElementById('resultText').innerHTML = "<p>There was a server error. Please try again later.</p>";
    } else { // Other Errors (catch-all)
        document.getElementById('resultText').innerHTML = "<p>An error occurred while processing the request. Please try again.</p>";
    }

    document.getElementById('warningContainer').style.display = 'none';
    document.getElementById('resultContainer').style.display = 'block';
    console.error("Error details (status code " + response.status + "):", response.statusText);
}

function toggleLanguage() {
    const currentLang = document.documentElement.lang === 'ar' ? 'en' : 'ar';
    document.documentElement.lang = currentLang;
    localStorage.setItem('preferredLanguage', currentLang);

    updateLanguage();
}

function updateLanguage() {
    const currentLang = document.documentElement.lang;

    document.querySelectorAll('[data-key]').forEach(function (element) {
        const key = element.getAttribute('data-key');
        element.textContent = translations[currentLang][key];
    });

    // Update placeholders
    document.querySelectorAll('input[placeholder]').forEach(function (input) {
        const key = input.getAttribute('data-key');
        input.setAttribute('placeholder', translations[currentLang][key]);
    });
}

document.addEventListener('DOMContentLoaded', function () {
    const preferredLanguage = localStorage.getItem('preferredLanguage') || 'en';
    document.documentElement.lang = preferredLanguage;
    updateLanguage();
});
