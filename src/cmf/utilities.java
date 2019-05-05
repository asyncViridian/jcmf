package cmf;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.TimeUnit;
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

    public static boolean fileExists(String filePathString) {
        boolean r = false;
        File f = new File(filePathString);
        if (f.exists() && !f.isDirectory()) {
            r = true;
        }
        return r;
    }

    //note filename should give full path
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
    public static synchronized cmdOut runCmd(String[] args, int timeout) {
        String ret = Arrays.stream(args).collect(joining(" ")) + "executing";
        boolean exitVal = true;
        ArrayList<String> out = new ArrayList<>();
        try {
            ProcessBuilder probuilder = new ProcessBuilder(args);
            Process p = probuilder.start();
            // output to java console
            BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String line;
            while ((line = reader.readLine()) != null) {
                out.add(line);
                System.out.println(line);
            }
            exitVal = p.waitFor(timeout, TimeUnit.SECONDS); //10 minute timeout guard. Java 8 only
            if (exitVal) {
                ret = ret + " command process completed successfully.";
            } else {
                ret = ret + " command process is killed over " + timeout + " seconds";
            }
        } catch (IOException | InterruptedException ex) {
            System.out.println(ex.getMessage());
        }
        System.out.println(ret);
        return new cmdOut(out, exitVal);
    }

    public static synchronized cmdOut runCmd(String[] args) {
        String ret = Arrays.stream(args).collect(joining(" ")) + "executing";
        int exitCode = 0;
        boolean res = true;
        ArrayList<String> out = new ArrayList<>();
        try {
            ProcessBuilder probuilder = new ProcessBuilder(args);
            Process p = probuilder.start();
            // output to java console
            BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String line;
            while ((line = reader.readLine()) != null) {
                out.add(line);
                System.out.println(line);
            }
            exitCode = p.waitFor(); //no time limit
            if (exitCode == 0) {
                ret = ret + " command process completed successfully.";
                res = true;
            } else {
                ret = ret + " command process failed.";
                // couldn't produce acceptable output
                res = false;
            }
        } catch (IOException | InterruptedException ex) {
            System.out.println(ex.getMessage());
        }
        System.out.println(ret);
        return new cmdOut(out, res, exitCode);
    }

}
