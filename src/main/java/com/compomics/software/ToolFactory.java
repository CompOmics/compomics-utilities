package com.compomics.software;

import com.compomics.software.cli.CommandLineUtils;
import com.compomics.software.dialogs.PeptideShakerSetupDialog;
import com.compomics.software.dialogs.ReporterSetupDialog;
import com.compomics.software.dialogs.SearchGuiSetupDialog;
import com.compomics.util.parameters.UtilitiesUserParameters;
import java.io.*;
import java.util.ArrayList;
import javax.swing.JFrame;
import javax.swing.JOptionPane;

/**
 * This class can be used to start compomics tools.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class ToolFactory {

    /**
     * The command line argument for a psdb file for PeptideShaker.
     */
    public static final String PEPTIDE_SHAKER_FILE_OPTION = "-psdb";
    /**
     * The command line argument for a zipped psdb URL for PeptideShaker.
     */
    public static final String PEPTIDE_SHAKER_URL_OPTION = "-zipUrl";
    /**
     * The command line argument for the download folder for the URL for
     * PeptideShaker.
     */
    public static final String PEPTIDE_SHAKER_URL_DOWNLOAD_FOLDER_OPTION = "-zipUrlFolder";
    /**
     * The command line argument to open a given PX accession in PRIDE Reshake.
     */
    public static final String PEPTIDE_SHAKER_PX_ACCESSION_OPTION = "-pxAccession";
    /**
     * The command line argument to indicate that the PX accession to open in
     * PRIDE Reshake is private. If not set public is assumed.
     */
    public static final String PEPTIDE_SHAKER_PX_ACCESSION_PRIVATE_OPTION = "-pxAccessionPrivate";
    /**
     * The command line argument for mgf files for SearchGUI.
     */
    public static final String SEARCHGUI_SPECTRUM_FILE_OPTION = "-mgf";
    public static final String searchGuiSpectrumFileOption = SEARCHGUI_SPECTRUM_FILE_OPTION;
    /**
     * The command line argument for FASTA file for SearchGUI.
     */
    public static final String SEARCHGUI_FASTA_FILE_OPTION = "-fasta";
    /**
     * The command line argument for raw files for SearchGUI.
     */
    public static final String SEARCHGUI_RAW_FILE_OPTION = "-raw";
    /**
     * The command line argument for a parameters file for SearchGUI.
     */
    public static final String SEARCHGUI_PARAMETERS_FILE_OPTION = "-identification_parameters";
    public static final String searchGuiParametersFileOption = SEARCHGUI_PARAMETERS_FILE_OPTION;
    /**
     * The command line argument for an output folder.
     */
    public static final String OUTPUT_FOLDER_OPTION = "-output_folder";
    public static final String outputFolderOption = OUTPUT_FOLDER_OPTION;
    /**
     * The command line argument for the species.
     */
    public static final String SPECIES_OPTION = "-species";
    /**
     * The command line argument for the species type.
     */
    public static final String SPECIES_TYPE_OPTION = "-species_type";
    /**
     * The command line argument for the PeptideShaker project name.
     */
    public static final String PROJECT_NAME_OPTION = "-project_name";
    /**
     * The command line argument for the config folder.
     */
    public static final String CONFIG_FOLDER = "-config_folder";

    /**
     * Starts PeptideShaker from the location of utilities preferences.
     *
     * @param parent a frame to display the path setting dialog (can be null)
     *
     * @throws IOException if an exception occurs while reading or writing a
     * file
     * @throws ClassNotFoundException if an exception occurs while reading the
     * user preferences
     * @throws InterruptedException if a threading issue occurs
     */
    public static void startPeptideShaker(JFrame parent) throws IOException, ClassNotFoundException, InterruptedException {
        startPeptideShaker(parent, null);
    }

    /**
     * Starts PeptideShaker from the location of utilities preferences and opens
     * the file given as argument. If null is given as file or if the file to
     * open is not found, the tool will go for a default start.
     *
     * @param parent a frame to display the path setting dialog (can be null)
     * @param psdbFile the file to open (psdb format) (can be null)
     *
     * @throws IOException if an exception occurs while reading or writing a
     * file
     * @throws ClassNotFoundException if an exception occurs while reading the
     * user preferences
     * @throws InterruptedException if a threading issue occurs
     */
    public static void startPeptideShaker(JFrame parent, File psdbFile) throws IOException, ClassNotFoundException, InterruptedException {

        UtilitiesUserParameters utilitiesUserPreferences = UtilitiesUserParameters.loadUserParameters();
        boolean openPeptideShaker = true;

        if (utilitiesUserPreferences.getPeptideShakerPath() == null || !(new File(utilitiesUserPreferences.getPeptideShakerPath()).exists())) {
            PeptideShakerSetupDialog peptideShakerSetupDialog = new PeptideShakerSetupDialog(parent, true);
            utilitiesUserPreferences = UtilitiesUserParameters.loadUserParameters();
            openPeptideShaker = !peptideShakerSetupDialog.isDialogCanceled();
        }
        if (openPeptideShaker) {
            if (utilitiesUserPreferences.getPeptideShakerPath() != null
                    && new File(utilitiesUserPreferences.getPeptideShakerPath()).exists()) {
                if (psdbFile != null) {
                    ArrayList<String> args = new ArrayList<>();
                    args.add(PEPTIDE_SHAKER_FILE_OPTION);
                    args.add(CommandLineUtils.getCommandLineArgument(psdbFile));
                    launch(utilitiesUserPreferences.getPeptideShakerPath(), "PeptideShaker", args);
                } else {
                    launch(utilitiesUserPreferences.getPeptideShakerPath(), "PeptideShaker");
                }
            } else {
                throw new IllegalArgumentException("PeptideShaker not found in " + utilitiesUserPreferences.getPeptideShakerPath());
            }
        }
    }

    /**
     * Starts PeptideShaker from the location of utilities preferences in the
     * Reshake mode and attempts at selecting the given project.
     *
     * @param parent a frame to display the path setting dialog (can be null)
     * @param pxAccession the ProteomeXchange accession of the project to open
     * (can be null)
     *
     * @throws IOException if an exception occurs while reading or writing a
     * file
     * @throws ClassNotFoundException if an exception occurs while reading the
     * user preferences
     * @throws InterruptedException if a threading issue occurs
     */
    public static void startReshake(JFrame parent, String pxAccession) throws IOException, ClassNotFoundException, InterruptedException {

        UtilitiesUserParameters utilitiesUserPreferences = UtilitiesUserParameters.loadUserParameters();
        boolean openPeptideShaker = true;

        if (utilitiesUserPreferences.getPeptideShakerPath() == null || !(new File(utilitiesUserPreferences.getPeptideShakerPath()).exists())) {
            PeptideShakerSetupDialog peptideShakerSetupDialog = new PeptideShakerSetupDialog(parent, true);
            utilitiesUserPreferences = UtilitiesUserParameters.loadUserParameters();
            openPeptideShaker = !peptideShakerSetupDialog.isDialogCanceled();
        }
        if (openPeptideShaker) {
            if (utilitiesUserPreferences.getPeptideShakerPath() != null
                    && new File(utilitiesUserPreferences.getPeptideShakerPath()).exists()) {
                if (pxAccession != null) {
                    ArrayList<String> args = new ArrayList<>();
                    args.add(PEPTIDE_SHAKER_PX_ACCESSION_OPTION);
                    args.add(pxAccession);
                    launch(utilitiesUserPreferences.getPeptideShakerPath(), "PeptideShaker", args);
                } else {
                    launch(utilitiesUserPreferences.getPeptideShakerPath(), "PeptideShaker");
                }
            } else {
                throw new IllegalArgumentException("PeptideShaker not found in " + utilitiesUserPreferences.getPeptideShakerPath());
            }
        }
    }

    /**
     * Starts PeptideShaker from the location of utilities preferences and opens
     * the file given as argument. If null is given as file or if the file to
     * open is not found, the tool will go for a default start.
     *
     * @param parent a frame to display the path setting dialog (can be null)
     * @param zipUrl the URL with the zipped PeptideShaker project to open (can
     * be null)
     * @param downloadUrlFolder the folder to download the project to, mandatory
     * if zipUrl is used
     *
     * @throws IOException if an exception occurs while reading or writing a
     * file
     * @throws ClassNotFoundException if an exception occurs while reading the
     * user preferences
     * @throws InterruptedException if a threading issue occurs
     */
    public static void startPeptideShakerFromURL(JFrame parent, String zipUrl, String downloadUrlFolder) throws IOException, ClassNotFoundException, InterruptedException {

        UtilitiesUserParameters utilitiesUserPreferences = UtilitiesUserParameters.loadUserParameters();
        boolean openPeptideShaker = true;

        if (utilitiesUserPreferences.getPeptideShakerPath() == null || !(new File(utilitiesUserPreferences.getPeptideShakerPath()).exists())) {
            PeptideShakerSetupDialog peptideShakerSetupDialog = new PeptideShakerSetupDialog(parent, true);
            utilitiesUserPreferences = UtilitiesUserParameters.loadUserParameters();
            openPeptideShaker = !peptideShakerSetupDialog.isDialogCanceled();
        }
        if (openPeptideShaker) {
            if (utilitiesUserPreferences.getPeptideShakerPath() != null
                    && new File(utilitiesUserPreferences.getPeptideShakerPath()).exists()) {
                if (zipUrl != null) {
                    ArrayList<String> args = new ArrayList<>();
                    args.add(PEPTIDE_SHAKER_URL_OPTION);
                    args.add(CommandLineUtils.getQuoteType() + zipUrl + CommandLineUtils.getQuoteType());
                    args.add(PEPTIDE_SHAKER_URL_DOWNLOAD_FOLDER_OPTION);
                    args.add(CommandLineUtils.getQuoteType() + downloadUrlFolder + CommandLineUtils.getQuoteType());
                    launch(utilitiesUserPreferences.getPeptideShakerPath(), "PeptideShaker", args);
                } else {
                    launch(utilitiesUserPreferences.getPeptideShakerPath(), "PeptideShaker");
                }
            } else {
                throw new IllegalArgumentException("PeptideShaker not found in " + utilitiesUserPreferences.getPeptideShakerPath());
            }
        }
    }

    /**
     * Starts Reporter from the location of utilities preferences.
     *
     * @param parent a frame to display the path setting dialog.
     *
     * @throws IOException if an exception occurs while reading or writing a
     * file
     * @throws ClassNotFoundException if an exception occurs while reading the
     * user preferences
     * @throws InterruptedException if a threading issue occurs
     */
    public static void startReporter(JFrame parent) throws IOException, ClassNotFoundException, InterruptedException {

        UtilitiesUserParameters utilitiesUserPreferences = UtilitiesUserParameters.loadUserParameters();
        boolean openReporter = true;

        if (utilitiesUserPreferences.getReporterPath() == null || !(new File(utilitiesUserPreferences.getReporterPath()).exists())) {
            ReporterSetupDialog reporterSetupDialog = new ReporterSetupDialog(parent, true);
            utilitiesUserPreferences = UtilitiesUserParameters.loadUserParameters();
            openReporter = !reporterSetupDialog.isDialogCanceled();
        }
        if (openReporter) {
            if (utilitiesUserPreferences.getReporterPath() != null
                    && new File(utilitiesUserPreferences.getReporterPath()).exists()) {
                launch(utilitiesUserPreferences.getReporterPath(), "Reporter");
            } else {
                throw new IllegalArgumentException("Reporter not found in " + utilitiesUserPreferences.getReporterPath());
            }
        }
    }

    /**
     * Starts SearchGUI from the location of utilities preferences.
     *
     * @param parent a frame to display the path setting dialog.
     *
     * @throws IOException if an exception occurs while reading or writing a
     * file
     * @throws ClassNotFoundException if an exception occurs while reading the
     * user preferences
     * @throws InterruptedException if a threading issue occurs
     */
    public static void startSearchGUI(JFrame parent) throws IOException, ClassNotFoundException, InterruptedException {
        startSearchGUI(parent, null, null, null, null, null, null, null, null, null);
    }

    /**
     * Starts SearchGUI from the location of utilities preferences.
     *
     * @param parent a frame to display the path setting dialog.
     * @param spectrumFiles the spectrum files to search (can be null)
     * @param rawFiles the raw files to search (can be null)
     * @param fastaFile the FASTA file
     * @param searchParameters the search parameters as a file (can be null)
     * @param outputFolder the output folder (can be null)
     * @param configFolder the config folder (can be null)
     * @param species the species (can be null)
     * @param speciesType the species type (can be null)
     * @param projectName the PeptideShaker project name
     *
     * @throws IOException if an exception occurs while reading or writing a
     * file
     * @throws ClassNotFoundException if an exception occurs while reading the
     * user preferences
     * @throws InterruptedException if a threading issue occurs
     */
    public static void startSearchGUI(
            JFrame parent,
            ArrayList<File> spectrumFiles,
            ArrayList<File> rawFiles,
            File fastaFile,
            File searchParameters,
            File outputFolder,
            File configFolder,
            String species,
            String speciesType,
            String projectName
    )
            throws IOException, ClassNotFoundException, InterruptedException {

        UtilitiesUserParameters utilitiesUserPreferences = UtilitiesUserParameters.loadUserParameters();
        boolean openSearchGUI = true;

        if (utilitiesUserPreferences.getSearchGuiPath() == null
                || !(new File(utilitiesUserPreferences.getSearchGuiPath()).exists())) {

            SearchGuiSetupDialog searchGuiSetupDialog = new SearchGuiSetupDialog(parent, true);
            utilitiesUserPreferences = UtilitiesUserParameters.loadUserParameters();
            openSearchGUI = !searchGuiSetupDialog.isDialogCanceled();

        }

        if (openSearchGUI) {

            if (utilitiesUserPreferences.getSearchGuiPath() != null
                    && new File(utilitiesUserPreferences.getSearchGuiPath()).exists()) {

                if (spectrumFiles == null && rawFiles == null && searchParameters == null && species == null) {

                    launch(utilitiesUserPreferences.getSearchGuiPath(), "SearchGUI");

                } else {

                    ArrayList<String> args = new ArrayList<>();

                    if (spectrumFiles != null && !spectrumFiles.isEmpty()) {
                        args.add(SEARCHGUI_SPECTRUM_FILE_OPTION);
                        args.add(CommandLineUtils.getCommandLineArgument(spectrumFiles));
                    }

                    if (rawFiles != null && !rawFiles.isEmpty()) {
                        args.add(SEARCHGUI_RAW_FILE_OPTION);
                        args.add(CommandLineUtils.getCommandLineArgument(rawFiles));
                    }

                    if (fastaFile != null) {
                        args.add(SEARCHGUI_FASTA_FILE_OPTION);
                        args.add(CommandLineUtils.getCommandLineArgument(fastaFile));
                    }

                    if (searchParameters != null) {
                        args.add(SEARCHGUI_PARAMETERS_FILE_OPTION);
                        args.add(CommandLineUtils.getCommandLineArgument(searchParameters));
                    }

                    if (outputFolder != null) {
                        args.add(OUTPUT_FOLDER_OPTION);
                        args.add(CommandLineUtils.getCommandLineArgument(outputFolder));
                    }

                    if (species != null) {
                        args.add(SPECIES_OPTION);
                        args.add(species);
                    }

                    if (speciesType != null) {
                        args.add(SPECIES_TYPE_OPTION);
                        args.add(speciesType);
                    }

                    if (projectName != null) {
                        args.add(PROJECT_NAME_OPTION);
                        args.add(projectName);
                    }

                    if (configFolder != null) {
                        args.add(CONFIG_FOLDER);
                        args.add(CommandLineUtils.getCommandLineArgument(configFolder));
                    }

                    launch(utilitiesUserPreferences.getSearchGuiPath(), "SearchGUI", args);

                }
            } else {
                throw new IllegalArgumentException("SearchGUI not found in " + utilitiesUserPreferences.getSearchGuiPath());
            }
        }
    }

    /**
     * Generic method to start a tool.
     *
     * @param toolPath the path to the tool
     * @param toolName the name of the tool (used for bug report)
     * @param args the arguments to pass to the tool (ignored if null)
     */
    private static void launch(String toolPath, String toolName) throws IOException, InterruptedException {
        launch(toolPath, toolName, null);
    }

    /**
     * Generic method to start a tool.
     *
     * @param toolPath the path to the tool
     * @param toolName the name of the tool (used for bug report)
     * @param args the arguments to pass to the tool (ignored if null)
     */
    private static void launch(String toolPath, String toolName, ArrayList<String> args) throws IOException, InterruptedException {

        String arguments = "";
        if (args != null) {
            arguments = CommandLineUtils.concatenate(args);
        }

        boolean debug = false;

        String quote = CommandLineUtils.getQuoteType();

        String cmdLine = "java -jar " + quote + toolPath + quote + " " + arguments;

        if (debug) {
            System.out.println(cmdLine);
        }

        Process p = Runtime.getRuntime().exec(cmdLine);

        InputStream stderr = p.getErrorStream();
        InputStreamReader isr = new InputStreamReader(stderr);
        BufferedReader br = new BufferedReader(isr);

        String temp = "<ERROR>" + System.getProperty("line.separator") + System.getProperty("line.separator");

        if (debug) {
            System.out.println("<ERROR>");
        }

        String line = br.readLine();

        boolean error = false;

        while (line != null) {

            if (debug) {
                System.out.println(line);
            }

            temp += line + System.getProperty("line.separator");
            line = br.readLine();
            error = true;
        }

        br.close();
        isr.close();
        stderr.close();

        if (debug) {
            System.out.println("</ERROR>");
        }

        temp += System.getProperty("line.separator") + "The command line executed:" + System.getProperty("line.separator");
        temp += cmdLine + System.getProperty("line.separator");
        temp += System.getProperty("line.separator") + "</ERROR>" + System.getProperty("line.separator");
        int exitVal = p.waitFor();

        if (debug) {
            System.out.println("Process exitValue: " + exitVal);
        }

        if (error) {
            File logFile = new File(new File(toolPath).getParentFile().getAbsolutePath() + "/resources/", toolName + ".log");
            FileWriter f = new FileWriter(logFile, true);
            f.write(System.getProperty("line.separator") + System.getProperty("line.separator") + temp + System.getProperty("line.separator") + System.getProperty("line.separator"));
            f.close();

            javax.swing.JOptionPane.showMessageDialog(null,
                    "Failed to start " + toolName + ".\n\n"
                    + "Inspect the log file for details: resources/" + toolName + ".log.\n\n"
                    + "Then go to Troubleshooting at https://compomics.github.io/projects/peptide-shaker.html.",
                    toolName + " - Startup Failed", JOptionPane.ERROR_MESSAGE);
        }
    }
}
