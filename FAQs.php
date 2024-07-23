<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pharmacare+ FAQ</title>
    <link rel="stylesheet" href="home2.css">
    <style>
        .faq-container {
            text-align: left;
            margin: 20px auto;
            max-width: 800px;
            background-color: #fff;
            padding: 20px;
            box-shadow: 0 0 10px rgba(0, 0, 0, 0.1);
            border-radius: 10px;
        }

        .faq-container h2 {
            color: #333;
            font-size: 24px;
            border-bottom: 1px solid #ddd;
            padding-bottom: 10px;
        }

        .faq-container h3 {
            color: #008080;
            font-size: 20px;
            margin-top: 20px;
        }

        .faq-container p {
            color: #555;
            font-size: 16px;
            line-height: 1.6;
        }

        .faq-item {
            margin-bottom: 20px;
        }

        .faq-item button {
            background: none;
            border: none;
            font-size: 18px;
            color: #008080;
            cursor: pointer;
            text-align: left;
            width: 100%;
            padding: 10px;
            border-bottom: 1px solid #ddd;
        }

        .faq-item button:focus {
            outline: none;
        }

        .faq-answer {
            display: none;
            padding: 10px 0;
        }

        .faq-item button.active + .faq-answer {
            display: block;
        }
    </style>
</head>
<body>
    <?php include "header.php" ?>
    <main>
        <button id="language-toggle">Toggle Language</button>
        <div class="faq-container">
            <h2 data-key="faq_title">Frequently Asked Questions</h2>
            
            <div class="faq-item">
                <button data-key="q1">What is the purpose of this service?</button>
                <div class="faq-answer">
                    <p data-key="a1">The purpose is to provide accurate information and warnings about drug interactions to help users avoid potential health problems.</p>
                </div>
            </div>

            <div class="faq-item">
                <button data-key="q2">How does the intelligent system provide warnings?</button>
                <div class="faq-answer">
                    <p data-key="a2">The system uses AI technologies to analyze drug databases and their potential interactions, then provides recommendations and warnings based on this analysis.</p>
                </div>
            </div>

            <div class="faq-item">
                <button data-key="q3">Can the intelligent system replace a doctor's consultation?</button>
                <div class="faq-answer">
                    <p data-key="a3">No, this service aims to provide additional information but does not replace consulting a doctor or pharmacist.</p>
                </div>
            </div>

            <div class="faq-item">
                <button data-key="q4">How can I use the service?</button>
                <div class="faq-answer">
                    <p data-key="a4">You can enter the names of the drugs you are taking or planning to take into the system, and it will analyze possible interactions between them and provide warnings and recommendations.</p>
                </div>
            </div>

            <div class="faq-item">
                <button data-key="q5">Is the drug database updated regularly?</button>
                <div class="faq-answer">
                    <p data-key="a5">Yes, the database is continuously updated to ensure the provided information is accurate and reliable.</p>
                </div>
            </div>

            <div class="faq-item">
                <button data-key="q6">Are my personal information and the drugs I take saved?</button>
                <div class="faq-answer">
                    <p data-key="a6">We respect your privacy. The information you enter is used only to analyze interactions and is not saved or shared with any third party.</p>
                </div>
            </div>

            <div class="faq-item">
                <button data-key="q7">Is there a cost to use this service?</button>
                <div class="faq-answer">
                    <p data-key="a7">The cost depends on the subscription type you choose. We offer both free and paid services with additional features.</p>
                </div>
            </div>

            <div class="faq-item">
                <button data-key="q8">Can I get detailed reports on drug interactions?</button>
                <div class="faq-answer">
                    <p data-key="a8">Yes, the system provides detailed reports that can be printed or shared with your doctor.</p>
                </div>
            </div>

            <div class="faq-item">
                <button data-key="q9">Does the service support languages other than Arabic?</button>
                <div class="faq-answer">
                    <p data-key="a9">Currently, the service is available in Arabic, but we are working on providing it in other languages soon.</p>
                </div>
            </div>

            <div class="faq-item">
                <button data-key="q10">How can I contact the support team?</button>
                <div class="faq-answer">
                    <p data-key="a10">You can contact us via email or phone listed on the contact page.</p>
                </div>
            </div>
        </div>
    </main>
    <?php include "footer.php" ?>
    <script src="translations.js"></script>
    <script src="main.js"></script>
</body>
</html>
