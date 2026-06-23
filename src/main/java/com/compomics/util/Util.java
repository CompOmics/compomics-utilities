package com.compomics.util;

import com.compomics.util.gui.waiting.waitinghandlers.ProgressDialogX;
import com.compomics.util.gui.file_handling.FileChooserUtil;
import java.awt.Color;
import java.awt.Component;
import java.io.*;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import javax.swing.JTable;

/**
 * Includes general help methods that are used by the other classes.
 *
 * @author Harald Barsnes
 * @author Marc Vaudel
 */
public class Util {

    /**
     * Forbidden characters in file names.
     */
    public static final String[] FORBIDDEN_CHARACTERS = {"!", ":", ";", "\\?", "/", "\\\\", "\\*", "<", ">", "\"", "'", "\\|"};
    /**
     * Default encoding, cf the second rule.
     */
    public static final String ENCODING = "UTF-8";
    /**
     * Default column separator for text files.
     */
    public static final String DEFAULT_COLUMN_SEPARATOR = "\t";
    /**
     * The line separator.
     */
    public static final String LINE_SEPARATOR = System.getProperty("line.separator");

    /**
     * The mass added per amino acid as part of the reference mass when
     * converting a tolerance in Dalton to ppm.
     */
    public static final double MASS_PER_AA = 100.0;

    /**
     * Removes characters from a string.
     *
     * @param string the string of interest
     * @param subString the sub-string to remove
     *
     * @return a version without forbidden characters
     */
    public static String removeSubString(
            String string,
            String subString
    ) {

        String result;
        String[] split = string.split(subString);

        if (split.length > 1) {

            StringBuilder stringBuilder = new StringBuilder(string.length());

            for (String splitPart : split) {

                stringBuilder.append(splitPart);

            }

            result = stringBuilder.toString();

        } else {

            result = string;

        }

        return result;

    }

    /**
     * Function for sending event requests to Google analytics.
     *
     * @param ua tracking ID
     * @param action the action
     * @param label the label
     *
     * @return true if the update was successful
     */
    public static boolean sendGAUpdate(
            String ua,
            String action,
            String label
    ) {

        // the plain java way
        boolean returnVal = true;
        String collect_url = "https://www.google-analytics.com/collect";
        String post = "v=1&tid=" + ua + "&cid=35119a79-1a05-49d7-b876-bb88420f825b&uid=asuueffeqqss&t=event&ec=usage&ea=" + action + "&el=" + label;

        try {

            HttpURLConnection connection = (HttpURLConnection) new URL(collect_url).openConnection();
            connection.setRequestMethod("POST");
            connection.setConnectTimeout(3000);
            connection.setReadTimeout(1000);
            connection.setDoOutput(true);
            DataOutputStream wr = new DataOutputStream(connection.getOutputStream());
            wr.writeBytes(post);
            connection.getResponseCode();

        } catch (IOException ex) {

            returnVal = false;

        }

        return returnVal;

    }

    /**
     * Retrieves the version number set in the pom file.
     *
     * @return the version number of compomics-utilities
     */
    public static String getVersion() {

        java.util.Properties p = new java.util.Properties();

        try {

            InputStream is = Util.class.getClassLoader().getResourceAsStream("compomics-utilities.properties");
            p.load(is);

        } catch (IOException e) {

            e.printStackTrace();

        }

        return p.getProperty("compomics-utilities.version");
    }

    /**
     * Removes characters from a string.
     *
     * @param string the string of interest
     * @param forbiddenCharacter the character to remove
     *
     * @return a version without forbidden characters
     */
    public static String removeCharacters(
            String string,
            char forbiddenCharacter
    ) {

        StringBuilder sb = new StringBuilder(string.length());

        char[] stringChars = string.toCharArray();

        for (int i = 0; i < stringChars.length; i++) {

            char charAtI = stringChars[i];

            if (charAtI == forbiddenCharacter) {

                continue;

            }

            sb.append(charAtI);

        }

        return sb.toString();
    }

    /**
     * Removes characters from a string.
     *
     * @param string the string of interest
     * @param forbiddenCharacters the characters to remove
     *
     * @return a version without forbidden characters
     */
    public static String removeCharacters(
            String string,
            String[] forbiddenCharacters
    ) {

        String result = string;

        for (String fc : forbiddenCharacters) {

            String[] split = result.split(fc);

            if (split.length > 1) {

                result = Arrays.stream(split)
                        .collect(Collectors.joining());

            }
        }

        return result;
    }

    /**
     * Removes the forbidden characters from a string.
     *
     * @param string the string of interest
     * @return a version without forbidden characters
     */
    public static String removeForbiddenCharacters(
            String string
    ) {

        return removeCharacters(string, FORBIDDEN_CHARACTERS);

    }

    /**
     * Indicates whether a string contains characters forbidden in file names.
     *
     * @param string the string of interest
     * @return a boolean indicating whether a string contains characters
     * forbidden in file names
     */
    public static boolean containsForbiddenCharacter(
            String string
    ) {

        for (String forbiddenCharacter : FORBIDDEN_CHARACTERS) {

            if (string.contains(forbiddenCharacter)) {

                return true;

            }
        }

        return false;

    }

    /**
     * Replaces all characters equals to a by b in the given string.
     *
     * @param string the string
     * @param a a
     * @param b b
     *
     * @return a string equal to string with a replaced by b
     */
    public static String replaceAll(
            String string,
            char a,
            char b
    ) {

        char[] stringArray = string.toCharArray();

        for (int i = 0; i < stringArray.length; i++) {

            if (stringArray[i] == a) {

                stringArray[i] = b;

            }
        }

        return new String(stringArray);

    }

    /**
     * Rounds a double value to the wanted number of decimal places.
     *
     * @param d the double to round of
     * @param places number of decimal places wanted
     * @return double - the new double
     */
    public static double roundDouble(
            double d,
            int places
    ) {

        return Math.round(d * Math.pow(10, (double) places)) / Math.pow(10, (double) places);

    }

    /**
     * Floors a double value to the wanted number of decimal places.
     *
     * @param d the double to floor of
     * @param places number of decimal places wanted
     *
     * @return double - the new double
     */
    public static double floorDouble(
            double d,
            int places
    ) {

        return (new BigDecimal(String.valueOf(d)).setScale(places, RoundingMode.FLOOR)).doubleValue();

    }

    /**
     * Returns the ppm value of the given mass error relative to its theoretical
     * m/z value.
     *
     * @param theoreticalMzValue the theoretical mass
     * @param massError the mass error
     *
     * @return the mass error as a ppm value relative to the theoretical mass
     */
    public static double getPpmError(
            double theoreticalMzValue,
            double massError
    ) {

        double ppmValue = (massError / theoreticalMzValue) * 1000000;
        return ppmValue;

    }

    /**
     * Converts a color to hex format for use in HTML tags.
     *
     * @param colorRGB the color in RGB representation
     *
     * @return the color in hex format
     */
    public static String color2Hex(
            int colorRGB
    ) {

        return Integer.toHexString(colorRGB & 0x00ffffff);

    }

    /**
     * Converts a color to hex format for use in HTML tags.
     *
     * @param color the color to convert
     *
     * @return the color in hex format
     */
    public static String color2Hex(
            Color color
    ) {

        return color2Hex(color.getRGB());

    }

    /**
     * Returns the color object corresponding to the given RGB representation.
     *
     * @param colorRGB the color in RGBs representation
     *
     * @return the color object
     */
    public static Color getColor(
            int colorRGB
    ) {

        return new Color((colorRGB >> 16) & 0xFF, (colorRGB >> 8) & 0xFF, colorRGB & 0xFF);

    }

    /**
     * Returns the table as a separated text file.
     *
     * @param table the table to turn in to text
     * @param separator the text separator
     * @param progressDialog the progress dialog
     * @param removeHtml if true, HTML is converted to text
     *
     * @return the table as a separated text file
     */
    public static String tableToText(
            JTable table,
            String separator,
            ProgressDialogX progressDialog,
            boolean removeHtml
    ) {

        StringBuilder tableAsString = new StringBuilder();

        for (int i = 0; i < table.getColumnCount() && !progressDialog.isRunCanceled(); i++) {
            tableAsString.append(table.getColumnName(i)).append(separator);
        }

        progressDialog.setPrimaryProgressCounterIndeterminate(false);
        progressDialog.setMaxPrimaryProgressCounter(table.getRowCount());

        String lineBreak = System.getProperty("line.separator");

        tableAsString.append(lineBreak);

        for (int i = 0; i < table.getRowCount() && !progressDialog.isRunCanceled(); i++) {

            progressDialog.increasePrimaryProgressCounter();

            for (int j = 0; j < table.getColumnCount() && !progressDialog.isRunCanceled(); j++) {

                if (table.getValueAt(i, j) != null) {

                    String tempValue = table.getValueAt(i, j).toString();

                    // remove html tags
                    if (tempValue.contains("<html>") && removeHtml) {
                        tempValue = tempValue.replaceAll("\\<[^>]*>", "");
                    }

                    tableAsString.append(tempValue).append(separator);

                } else {

                    tableAsString.append(separator);

                }
            }

            tableAsString.append(lineBreak);

        }

        return tableAsString.toString();

    }

    /**
     * Writes the table to a file as separated text.
     *
     * @param table the table to write to file
     * @param separator the text separator
     * @param progressDialog the progress dialog
     * @param removeHtml if true, HTML is converted to text
     * @param writer the writer where the file is to be written
     *
     * @throws IOException if a problem occurs when writing to the file
     */
    public static void tableToFile(
            JTable table,
            String separator,
            ProgressDialogX progressDialog,
            boolean removeHtml,
            BufferedWriter writer
    ) throws IOException {

        for (int i = 0; i < table.getColumnCount(); i++) {

            writer.write(table.getColumnName(i) + separator);

        }

        if (progressDialog != null) {

            progressDialog.setPrimaryProgressCounterIndeterminate(false);
            progressDialog.setMaxPrimaryProgressCounter(table.getRowCount());

        }

        String lineBreak = System.getProperty("line.separator");

        writer.write(lineBreak);

        for (int i = 0; i < table.getRowCount(); i++) {

            if (progressDialog != null) {

                if (progressDialog.isRunCanceled()) {

                    return;

                }

                progressDialog.increasePrimaryProgressCounter();
            }

            for (int j = 0; j < table.getColumnCount(); j++) {

                if (progressDialog != null) {

                    if (progressDialog.isRunCanceled()) {

                        return;

                    }

                }

                if (table.getValueAt(i, j) != null) {

                    String tempValue = table.getValueAt(i, j).toString();

                    // remove html tags
                    if (tempValue.contains("<html>") && removeHtml) {

                        tempValue = tempValue.replaceAll("\\<[^>]*>", "");

                    }

                    writer.write(tempValue + separator);

                } else {

                    writer.write(separator);

                }

            }

            writer.write(lineBreak);
        }
    }

    /**
     * Convenience methods indicating whether the content of two lists have the
     * same content. Equality is based on the hash of the objects. Note that
     * this method does not work for lists containing null.
     *
     * @param list1 the first list
     * @param list2 the second list
     *
     * @return a boolean indicating whether list1 has the same content as list2
     */
    public static boolean sameLists(
            ArrayList<?> list1,
            ArrayList<?> list2
    ) {

        if (list1.size() != list2.size()) {

            return false;

        }

        HashMap<Object, Long> list1Occurrence = list1.stream()
                .collect(
                        Collectors.groupingBy(
                                a -> a,
                                HashMap::new,
                                Collectors.counting()
                        )
                );

        HashMap<Object, Long> list2Occurrence = list2.stream()
                .collect(
                        Collectors.groupingBy(
                                a -> a,
                                HashMap::new,
                                Collectors.counting()
                        )
                );

        return list1Occurrence.entrySet().stream()
                .allMatch(
                        entry -> list2Occurrence.containsKey(entry.getKey()) && list2Occurrence.get(entry.getKey()).equals(entry.getValue())
                );

    }

    /**
     * Returns the occurrence of a character in a string.
     *
     * @param input the string of interest
     * @param character the character to look for
     *
     * @return the occurrence of a character in a string
     */
    public static int getOccurrence(
            String input,
            char character
    ) {

        return (int) input.chars()
                .filter(aa -> aa == character)
                .count();

    }

    /**
     * Returns at which indexes a small string can be found in a big string.
     *
     * @param bigString the big string
     * @param smallString the small string
     * @return a list of the indexes where the small string can be found in the
     * big string
     */
    public static ArrayList<Integer> getIndexes(
            String bigString,
            String smallString
    ) {
        
        Pattern pattern = Pattern.compile(smallString);
        ArrayList<Integer> result = new ArrayList<>();
        Matcher matcher = pattern.matcher(bigString);
        matcher.matches();
        int index = 0;

        while (matcher.find(index)) {

            index = matcher.start();
            index++;
            result.add(index);

        }

        return result;

    }

    /**
     * Method for reading a double value as a string which uses either "," or
     * "." as the decimal symbol.
     *
     * @param doubleAsString the double value as a string
     * @return the double value
     * @throws NumberFormatException thrown if the double cannot be read as a
     * double
     */
    public static double readDoubleAsString(
            String doubleAsString
    ) throws NumberFormatException {

        BigDecimal temp;

        try {

            temp = new BigDecimal(doubleAsString);

        } catch (NumberFormatException e) {

            doubleAsString = doubleAsString.replaceAll("\\.", "");
            doubleAsString = doubleAsString.replaceAll(",", "\\.");

            try {

                temp = new BigDecimal(doubleAsString);

            } catch (NumberFormatException ex) {

                throw new NumberFormatException(doubleAsString + " cannot be read as a floating value!");

            }
        }

        return temp.doubleValue();

    }

    /**
     * Converts a boolean value to the corresponding integer value, 0 for false
     * and 1 for true.
     *
     * @param booleanToConvert the boolean value to convert
     *
     * @return 0 for false and 1 for true
     */
    public static int convertBooleanToInteger(
            Boolean booleanToConvert
    ) {

        return booleanToConvert ? 1 : 0;

    }

    /**
     * Returns a string in the form key(value).
     *
     * @param key the key
     * @param value the value
     *
     * @return a string in the form value(attribute)
     */
    public static String keyValueToString(
            String key,
            String value
    ) {

        StringBuilder sb = new StringBuilder(key.length() + value.length() + 2);
        sb.append(key).append("(").append(value).append(")");

        return sb.toString();

    }

    public static File getUserSelectedFile(
            Component parent,
            String fileEnding,
            String fileFormatDescription,
            String dialogTitle,
            String lastSelectedFolder,
            String suggestedFileName,
            boolean openDialog
    ) {
        return FileChooserUtil.getUserSelectedFile(
                parent,
                fileEnding,
                fileFormatDescription,
                dialogTitle,
                lastSelectedFolder,
                suggestedFileName,
                openDialog
        );
    }

    public static boolean deleteDir(File directory) {
        if (directory == null || !directory.exists()) {
            return true;
        }
        File[] files = directory.listFiles();
        if (files != null) {
            for (File file : files) {
                if (file.isDirectory()) {
                    deleteDir(file);
                } else if (!file.delete()) {
                    return false;
                }
            }
        }
        return directory.delete();
    }

    public static String getFileName(File file) {
        String name = file.getName();
        int extensionIndex = name.lastIndexOf('.');
        return extensionIndex > 0 ? name.substring(0, extensionIndex) : name;
    }
}
