package cmf;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.util.stream.Collectors.joining;

/**
 *
 * Class to execute external command
 */
public class utilities {

    // this is java 7 method, not functional,  also Java 8 Files.walk() can use lamdba / stream
    public static ArrayList<Path> findFiles(String pathStr, String filePattern) {
        Path pathSource = Paths.get(pathStr);
        MyFileVisitor finder = new MyFileVisitor(filePattern);
        try {
            Files.walkFileTree(pathSource, EnumSet.of(FileVisitOption.FOLLOW_LINKS), Integer.MAX_VALUE, finder);
        } catch (IOException ex) {
            System.err.println(ex);
        }
        return finder.getList();
    }

    //java 8 File.lis() is easier, not go to sub folders, it finds multiple files using filter startswith() foreach()
    public static String findFile(String path, String filename) {
        String f = null;
        try {
            f = Files.list(new File(path).toPath())
                    .filter(e -> e.getFileName().toString().equals(filename)).findAny().get().toString();
        } catch (IOException ex) {
            System.err.println(ex);
        }
        return f;
    }

    public static void deleteFile(String filename) {
        Path fileToDeletePath = Paths.get(filename);
        try {
            Files.delete(fileToDeletePath);
            System.out.println(filename + "deleted.");
        } catch (IOException ex) {
            System.err.println(ex);
        }
    }

    //return value is good to decide if the next command should be executed or not
    public static synchronized String runCmd(String[] args, int timeout) {
        String ret = Arrays.stream(args).collect(joining(" ")) + "executing";
        try {
            ProcessBuilder probuilder = new ProcessBuilder(args);
            Process p = probuilder.start();
            // output to java console
            BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String line;
            while ((line = reader.readLine()) != null) {
                System.out.println(line);
            }
            boolean exitVal = p.waitFor(timeout, TimeUnit.SECONDS); //10 minute timeout guard. 
            if (exitVal) {
                ret = Arrays.stream(args).collect(joining(" "))
                        + " command process completed successfully.";
            } else {
                ret = Arrays.stream(args).collect(joining(" "))
                        + " command process is killed over " + timeout + " seconds";
            }
        } catch (IOException | InterruptedException ex) {
            System.out.println(ex.getMessage());
        }
        return ret;
    }

}
