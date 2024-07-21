<?php include "header.php" ?>
    <main>
        <section class="interaction-checker">
            <h2 data-key="interaction_checker_title">Interaction Checker</h2>
            <div class="checker-container">
                <div class="form-container">
                    <h3 data-key="check_potential_risks">Check for potential risks of drug interaction :</h3>
                    <form id="interactionForm">
                        <div style="position: relative;">
                            <input type="text" placeholder="commercial name..." id="drug1" onkeyup="fetchSuggestions(this, 'suggestions1')" data-key="drug1_placeholder">
                            <div id="suggestions1" class="suggestions"></div>
                        </div>
                        <div style="position: relative;">
                            <input type="text" placeholder="commercial name..." id="drug2" onkeyup="fetchSuggestions(this, 'suggestions2')" data-key="drug2_placeholder">
                            <div id="suggestions2" class="suggestions"></div>
                        </div>
                        <p data-key="add_two_drugs">Add two drugs.</p>
                        <div class="button-group">
                            <button type="button" onclick="checkInteraction()" data-key="check_interaction">Check interaction</button>
                            <button type="reset" onclick="clearResult()" data-key="clear">Clear</button>
                        </div>
                    </form>
                </div>
                <div class="warning-container" id="warningContainer">
                    <h3><span class="warning-icon">⚠️</span><span data-key="warning_title"> Warning :</span></h3>
                    <p data-key="warning_message">The results of checking are based on the current knowledge. If no interactions are found between two drugs, it does not necessarily mean that no interactions exist. Always consult with a healthcare professional.</p>
                </div>
                <div class="result-container" id="resultContainer" style="display: none;">
                    <h3 data-key="interaction_result">Interaction Result:</h3>
                    <p id="resultText"></p>
                </div>
            </div>
        </section>
    </main>
    <?php include "footer.php" ?>