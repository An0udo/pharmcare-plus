<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pharmacare+ Explanation</title>
    <link rel="stylesheet" href="home2.css">
    <style>
        .explanation-container {
            text-align: left;
            margin: 20px auto;
            max-width: 800px;
            background-color: #fff;
            padding: 20px;
            box-shadow: 0 0 10px rgba(0, 0, 0, 0.1);
            border-radius: 10px;
        }

        .explanation-container h2 {
            color: #333;
            font-size: 24px;
            border-bottom: 1px solid #ddd;
            padding-bottom: 10px;
        }

        .explanation-container h3 {
            color: #008080;
            font-size: 20px;
            margin-top: 20px;
        }

        .explanation-container p {
            color: #555;
            font-size: 16px;
            line-height: 1.6;
        }

        .icon-container {
            display: flex;
            justify-content: space-between;
            margin-top: 20px;
        }

        .icon-container div {
            text-align: center;
            flex: 1;
        }

        .icon-container h4 {
            color: #333;
            font-size: 16px;
        }

        .responsive-img {
            max-width: 100%;
            height: auto;
            margin: 20px 0;
        }
    </style>
</head>
<body>
    <?php include "header.php" ?>
    <main>
        <div class="explanation-container">
            <h2 data-key="explanation_title">Explanation</h2>
            <h3 data-key="what_is_pharmacare">What is Pharmacare+?</h3>
            <p data-key="pharmacare_description">Pharmacare+ relies on sophisticated techniques in the fields of analytical analysis and artificial intelligence to enhance drug safety and patient care. The application analyzes drug interaction data using machine learning models, enabling it to accurately identify potential interactions between drugs and conflicts.</p>
            
            <h3 data-key="what_can_ddinter_do">What can DDInter do?</h3>
            <p data-key="ddinter_description">This application will serve as a comprehensive database system designed to utilize Artificial Intelligence for the Detection and Management of Drug Interactions, enhancing patient safety and healthcare efficiency in the country.</p>
            
            <img src="Screenshot 2024-07-06 at 4.35.00 PM.png" alt="Explanation" class="responsive-img">
            
            <div class="icon-container">
                <div>
                    <h4 data-key="interaction_description">Interaction description</h4>
                </div>
                <div>
                    <h4 data-key="severity_level">Severity level</h4>
                </div>
                <div>
                    <h4 data-key="management_methods">Management methods</h4>
                </div>
                <div>
                    <h4 data-key="references">References</h4>
                </div>
            </div>
        </div>
    </main>
    <?php include "footer.php" ?>
    