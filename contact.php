<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pharmacare+ Contact Us</title>
    <link rel="stylesheet" href="home2.css">
    <style>
        .contact-container {
            text-align: left;
            margin: 20px auto;
            max-width: 800px;
            padding: 20px;
        }

        .contact-container h2 {
            color: #333;
            font-size: 24px;
            border-bottom: 1px solid #ddd;
            padding-bottom: 10px;
        }

        .contact-container p {
            color: #555;
            font-size: 16px;
            line-height: 1.6;
        }

        .contact-form {
            display: flex;
            flex-direction: column;
        }

        .contact-form input, .contact-form textarea {
            margin-bottom: 10px;
            padding: 10px;
            font-size: 16px;
            border: 1px solid #ddd;
            border-radius: 5px;
        }

        .contact-form button {
            padding: 10px;
            font-size: 18px;
            color: #fff;
            background-color: #333;
            border: none;
            border-radius: 5px;
            cursor: pointer;
        }

        .contact-form button:hover {
            background-color: #555;
        }
    </style>
</head>
<body>
    <?php include "header.php" ?>

    <main>
        <div class="contact-container">
            <h2>Let's Talk About Everything!</h2>
            <p>Lorem ipsum dolor sit amet, consectetur adipiscing elit. Voluptas debitis, fugit natus?</p>

            <?php
            if ($_SERVER["REQUEST_METHOD"] == "POST") {
                $name = htmlspecialchars(trim($_POST["name"]));
                $email = htmlspecialchars(trim($_POST["email"]));
                $subject = htmlspecialchars(trim($_POST["subject"]));
                $message = htmlspecialchars(trim($_POST["message"]));

                $to = "sannody0@outlook.sa"; // Replace with the email address you want to send messages to
                $headers = "From: $email\r\n";
                $headers .= "Reply-To: $email\r\n";
                $headers .= "Content-Type: text/html; charset=UTF-8\r\n";

                $full_message = "<html><body>";
                $full_message .= "<h2>New message from Pharmacare+ website</h2>";
                $full_message .= "<p><strong>Name:</strong> $name</p>";
                $full_message .= "<p><strong>Email:</strong> $email</p>";
                $full_message .= "<p><strong>Subject:</strong> $subject</p>";
                $full_message .= "<p><strong>Message:</strong><br>" . nl2br($message) . "</p>";
                $full_message .= "</body></html>";

                if (mail($to, $subject, $full_message, $headers)) {
                    echo "<p>Your message has been sent successfully.</p>";
                } else {
                    echo "<p>Sorry, there was an error sending your message.</p>";
                }
            }
            ?>

            <form class="contact-form" method="POST" action="">
                <input type="text" name="name" placeholder="Name" required>
                <input type="email" name="email" placeholder="Email" required>
                <input type="text" name="subject" placeholder="Subject" required>
                <textarea name="message" rows="5" placeholder="Leave your message" required></textarea>
                <button type="submit">Send</button>
            </form>
        </div>
    </main>
    
    <?php include "footer.php" ?>
</body>
</html>
