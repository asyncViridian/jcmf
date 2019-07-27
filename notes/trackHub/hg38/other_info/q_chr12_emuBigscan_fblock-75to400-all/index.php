<?php
echo "<html><body style=\"font-family:monospace;\">";
$path = ".";
$dh = opendir($path);
$i=1;
$files=scandir($path);
foreach ($files as $file) {
    if($file != "." && $file != ".." && $file != "index.php" && $file != ".htaccess" && $file != "error_log" && $file != "cgi-bin") {
        echo "<a href='$path/$file'>$file</a><br /><br />";
        $i++;
    }
}
closedir($dh);
echo "</body></html>";
?> 