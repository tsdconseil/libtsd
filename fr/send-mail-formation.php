<?php
/**
* Author: Luis Zuno
* Email: luis@luiszuno.com
* URL: http://www.luiszuno.com
* Version: 1.0.0 
**/

//vars
$subject = $_POST['subject'];
$to = explode(',', $_POST['to'] );

$from = $_POST['email'];

//data
$msg = "NOM : "  .$_POST['name']    ."<br>\n";
$msg .= "COURRIEL : "  .$_POST['email']    ."<br>\n";
$msg .= "DOMAINES DE FORMATION : "  .$_POST['dom']    ."<br>\n";
$msg .= "SOCIETE : "  .$_POST['soc']    ."<br>\n";
$msg .= "SITE WEB : "  .$_POST['web']    ."<br>\n";
$msg .= "ADRESSE : "  .$_POST['adr']    ."<br>\n";
$msg .= "COMMENTAIRES : "  .$_POST['comments']    ."<br>\n";

//Headers
$headers  = "MIME-Version: 1.0\r\n";
$headers .= "Content-type: text/html; charset=UTF-8\r\n";
$headers .= "From: <".$from. ">" ;


//send for each mail
foreach($to as $mail){
   mail($mail, $subject, $msg, $headers);
}

?>
