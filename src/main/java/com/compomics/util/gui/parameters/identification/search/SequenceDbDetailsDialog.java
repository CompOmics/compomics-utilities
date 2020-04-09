package com.compomics.util.gui.parameters.identification.search;

import com.compomics.util.Util;
import com.compomics.util.examples.BareBonesBrowserLaunch;
import com.compomics.util.experiment.biology.proteins.Protein;
import com.compomics.util.experiment.biology.taxonomy.SpeciesFactory;
import com.compomics.util.experiment.io.biology.protein.FastaParameters;
import com.compomics.util.experiment.io.biology.protein.FastaSummary;
import com.compomics.util.experiment.io.biology.protein.Header;
import com.compomics.util.experiment.io.biology.protein.ProteinDatabase;
import com.compomics.util.experiment.io.biology.protein.converters.DecoyConverter;
import com.compomics.util.experiment.io.biology.protein.iterators.FastaIterator;
import com.compomics.util.gui.JOptionEditorPane;
import com.compomics.util.gui.protein.FastaParametersDialog;
import com.compomics.util.gui.waiting.waitinghandlers.ProgressDialogX;
import com.compomics.util.io.IoUtil;
import com.compomics.util.io.file.LastSelectedFolder;
import com.compomics.util.parameters.UtilitiesUserParameters;
import java.awt.Dialog;
import java.awt.Frame;
import java.awt.Image;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.SpinnerListModel;
import javax.swing.filechooser.FileFilter;

/**
 * This dialog displays information about a sequence database.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class SequenceDbDetailsDialog extends javax.swing.JDialog {

    /**
     * A simple progress dialog.
     */
    private static ProgressDialogX progressDialog;
    /**
     * The last selected folder.
     */
    private LastSelectedFolder lastSelectedFolder = null;
    /**
     * Boolean indicating whether the database can be changed.
     */
    private boolean dbEditable = true;
    /**
     * The icon to display when waiting.
     */
    private Image waitingImage;
    /**
     * The normal icon.
     */
    private Image normalImange;
    /**
     * The parent frame.
     */
    private Frame parentFrame;
    /**
     * The utilities user parameters.
     */
    private UtilitiesUserParameters utilitiesUserParameters;
    /**
     * The key to use to store FASTA files paths.
     */
    public static final String lastFolderKey = "fastaFile";
    /**
     * The selected FASTA file.
     */
    private String selectedFastaFile = null;
    /**
     * The parameters used to parse the FASTA file.
     */
    private FastaParameters fastaParameters = null;
    /**
     * Summary information on the FASTA file content.
     */
    private FastaSummary fastaSummary = null;
    /**
     * The batch size of proteins to sample.
     */
    public static final int SAMPLE_BATCH_SIZE = 100;
    /**
     * Accessions of the sampled proteins.
     */
    private ArrayList<String> accessionsSample = new ArrayList<>(SAMPLE_BATCH_SIZE);
    /**
     * Sample of proteins from the database.
     */
    private HashMap<String, Protein> proteinsSample = new HashMap<>(SAMPLE_BATCH_SIZE);
    /**
     * Sample of protein headers from the database.
     */
    private HashMap<String, Header> headersSample = new HashMap<>(SAMPLE_BATCH_SIZE);
    /**
     * A protein iterator to fill the sample.
     */
    private FastaIterator proteinIterator;
    /**
     * Boolean indicating whether the database selection was canceled by the
     * user.
     */
    private boolean canceled = false;

    /**
     * Creates a new SequenceDbDetailsDialog with a dialog as owner.
     *
     * @param owner the dialog owner
     * @param parent the parent frame
     * @param selectedFastaFile the selected FASTA file
     * @param fastaParameters the parameters used to parse the FASTA file
     * @param lastSelectedFolder the last selected folder
     * @param dbEditable if the database is editable
     * @param normalImange the normal icon
     * @param waitingImage the waiting icon
     */
    public SequenceDbDetailsDialog(Dialog owner, Frame parent, String selectedFastaFile, FastaParameters fastaParameters, LastSelectedFolder lastSelectedFolder, boolean dbEditable, Image normalImange, Image waitingImage) {

        super(owner, true);
        initComponents();
        this.parentFrame = parent;
        this.lastSelectedFolder = lastSelectedFolder;
        this.dbEditable = dbEditable;
        this.waitingImage = waitingImage;
        this.normalImange = normalImange;

        this.selectedFastaFile = selectedFastaFile;
        this.fastaParameters = fastaParameters;

        this.utilitiesUserParameters = UtilitiesUserParameters.loadUserParameters();

        editFastaParametersJLabel.setEnabled(dbEditable);

        if (this.selectedFastaFile != null) {

            loadFastaFile(selectedFastaFile, fastaParameters == null, dbEditable, true);

        }

        setLocationRelativeTo(owner);

    }

    /**
     * Creates a new SequenceDbDetailsDialog.
     *
     * @param parent the parent frame
     * @param selectedFastaFile the selected FASTA file
     * @param fastaParameters the parameters used to parse the FASTA file
     * @param lastSelectedFolder the last selected folder
     * @param dbEditable if the database is editable
     * @param normalImange the normal icon
     * @param waitingImage the waiting icon
     */
    public SequenceDbDetailsDialog(Frame parent, String selectedFastaFile, FastaParameters fastaParameters, LastSelectedFolder lastSelectedFolder, boolean dbEditable, Image normalImange, Image waitingImage) {

        super(parent, true);

        initComponents();

        this.parentFrame = parent;
        this.lastSelectedFolder = lastSelectedFolder;
        this.dbEditable = dbEditable;
        this.waitingImage = waitingImage;
        this.normalImange = normalImange;

        this.selectedFastaFile = selectedFastaFile;
        this.fastaParameters = fastaParameters;

        this.utilitiesUserParameters = UtilitiesUserParameters.loadUserParameters();

        editFastaParametersJLabel.setEnabled(dbEditable);

        if (this.selectedFastaFile != null) {

            loadFastaFile(selectedFastaFile, fastaParameters == null, dbEditable, true);

        }

        setLocationRelativeTo(parent);
    }

    /**
     * Set up the GUI.
     */
    private void setUpGUI() {

        if (selectedFastaFile != null) {

            File fastaFile = new File(selectedFastaFile);

            fileTxt.setText(selectedFastaFile);
            lastModifiedTxt.setText(new Date(fastaFile.lastModified()).toString());

            if (fastaParameters == null) {

                dbNameTxt.setText("");
                versionTxt.setText("");
                decoyButton.setEnabled(true);

            }

            if (fastaSummary != null) {

                dbNameTxt.setText(fastaSummary.getName());
                versionTxt.setText(fastaSummary.getVersion());

                decoyButton.setEnabled(!fastaSummary.containsDecoys() && dbEditable);

                // show the species present in the database
                speciesJTextField.setText(SpeciesFactory.getSpeciesDescription(fastaSummary.speciesOccurrence));

                // show the database type information
                HashMap<ProteinDatabase, Integer> databaseType = fastaSummary.databaseType;

                // the origin of the sequence information
                if (databaseType.size() == 1) {

                    ProteinDatabase proteinDatabase = databaseType.keySet().stream().findFirst().get();
                    typeJTextField.setText(proteinDatabase.getFullName());

                } else {

                    TreeMap<Integer, TreeSet<ProteinDatabase>> occurrenceToDBMap = databaseType.entrySet().stream()
                            .collect(Collectors.groupingBy(Entry::getValue)).entrySet().stream()
                            .collect(Collectors.toMap(
                                    Entry::getKey,
                                    entry -> entry.getValue().stream()
                                            .map(entry2 -> entry2.getKey())
                                            .collect(Collectors.toCollection(TreeSet::new)),
                                    (a, b) -> {
                                        a.addAll(b);
                                        return a;
                                    },
                                    TreeMap::new));

                    String dbOccurrenceText = occurrenceToDBMap.descendingMap().values().stream()
                            .flatMap(dbs -> dbs.stream())
                            .map(db -> db.getFullName() + " (" + databaseType.get(db) + ")")
                            .collect(Collectors.joining(", "));

                    typeJTextField.setText(dbOccurrenceText);

                }

                // the number of sequences
                String nSequences = fastaSummary.nSequences + " sequences";

                if (fastaSummary.containsDecoys()) {

                    nSequences += " (" + fastaSummary.nTarget + " target)";

                }

                sizeTxt.setText(nSequences);

            } else {

                speciesJTextField.setText("");
                typeJTextField.setText("");
                sizeTxt.setText("");

            }

            browseButton.setEnabled(dbEditable);

            if (fastaFile.exists()) {

                try {

                    proteinIterator = new FastaIterator(fastaFile);
                    bufferProteins();

                    accessionsSpinner.setEnabled(true);

                    updateSequence();

                } catch (Exception e) {

                    JOptionPane.showMessageDialog(this, "An error occurred while reading the fasta file.",
                            "Import error", JOptionPane.WARNING_MESSAGE);
                    e.printStackTrace();

                }

            } else {

                accessionsSpinner.setEnabled(false);

            }
        }
    }

    /**
     * Updates the displayed sequence.
     */
    private void updateSequence() {

        String accession = accessionsSpinner.getValue().toString();
        Header header = headersSample.get(accession);
        Protein protein = proteinsSample.get(accession);

        proteinTxt.setText(header.getRawHeader() + System.getProperty("line.separator") + protein.getSequence());
        proteinTxt.setCaretPosition(0);

    }

    /**
     * Returns the last selected folder.
     *
     * @return the last selected folder
     */
    public String getLastSelectedFolder() {
        if (lastSelectedFolder == null) {
            return null;
        }
        String folder = lastSelectedFolder.getLastSelectedFolder(lastFolderKey);
        if (folder == null) {
            folder = lastSelectedFolder.getLastSelectedFolder();
        }
        return folder;
    }

    /**
     * Allows the user to select a FASTA file, loads its information, and
     * returns a boolean indicating whether the process loading was successful.
     *
     * @param userCanDispose if true, the dialog is closed if the user cancels
     * the selection
     *
     * @return a boolean indicating whether a valid FASTA file was selected
     */
    public boolean selectDB(boolean userCanDispose) {

        File startLocation = null;

        if (utilitiesUserParameters.getDbFolder() != null && utilitiesUserParameters.getDbFolder().exists()) {

            startLocation = utilitiesUserParameters.getDbFolder();

        }

        if (startLocation == null) {

            startLocation = new File(getLastSelectedFolder());

        }

        JFileChooser fc = new JFileChooser(startLocation);

        FileFilter filter = new FileFilter() {

            @Override
            public boolean accept(File myFile) {

                return myFile.getName().toLowerCase().endsWith("fasta")
                        || myFile.isDirectory();
            }

            @Override
            public String getDescription() {
                return "FASTA (.fasta)";
            }

        };

        fc.setFileFilter(filter);
        int result = fc.showOpenDialog(this);

        if (result == JFileChooser.APPROVE_OPTION) {

            File file = fc.getSelectedFile();
            File folder = file.getParentFile();
            utilitiesUserParameters.setDbFolder(folder);
            lastSelectedFolder.setLastSelectedFolder(lastFolderKey, folder.getAbsolutePath());

            if (file.getName().contains(" ")) {

                file = renameFastaFileName(file);

                if (file == null) {

                    return false;

                }
            }

            loadFastaFile(file.getAbsolutePath(), fastaParameters == null, true, true);

            return true;

        } else if (userCanDispose) {

            dispose();
            return false;

        }

        return false;
    }

    /**
     * Loads the given FASTA file.
     *
     * @param fastaFile a FASTA file
     * @param inferParameters if true, FASTA parsing parameters are inferred
     * automatically
     * @param iNewFastaSummary if true, a new FASTA summary will be created even
     * if one already exists
     * @param setUpGUI if true the GUI will be updated
     */
    private void loadFastaFile(String fastaFile, boolean inferParameters, boolean iNewFastaSummary, boolean setUpGUI) {

        this.selectedFastaFile = fastaFile;
        final boolean newFastaSummary = iNewFastaSummary;

        progressDialog = new ProgressDialogX(this, parentFrame,
                normalImange,
                waitingImage,
                true);

        new Thread(new Runnable() {
            public void run() {
                try {
                    progressDialog.setVisible(true);
                } catch (IndexOutOfBoundsException e) {
                    // ignore
                }
            }
        }, "ProgressDialog").start();

        new Thread("importThread") {
            public void run() {

                try {

                    if (inferParameters) {

                        progressDialog.setPrimaryProgressCounterIndeterminate(true);
                        progressDialog.setTitle("Inferring Database Format. Please Wait...");
                        fastaParameters = FastaParameters.inferParameters(selectedFastaFile, progressDialog);

                    }

                    if (!progressDialog.isRunCanceled()) {

                        progressDialog.setWaitingText("Importing Database. Please Wait...");
                        fastaSummary = FastaSummary.getSummary(selectedFastaFile, fastaParameters, newFastaSummary, progressDialog);

                        if (!fastaSummary.containsDecoys()) {

                            int outcome = JOptionPane.showConfirmDialog(SequenceDbDetailsDialog.this,
                                    "The database does not seem to contain decoy sequences.\nAdd decoys?", "Add decoys?",
                                    JOptionPane.YES_NO_OPTION);

                            if (outcome == JOptionPane.YES_OPTION) {
                                generateTargetDecoyDatabase();
                            } else {
                                decoyButton.setEnabled(true);
                            }
                        }

                        if (setUpGUI && !progressDialog.isRunCanceled()) {
                            setUpGUI();
                        }

                        progressDialog.setRunFinished();

                    }

                } catch (IOException e) {
                    progressDialog.setRunFinished();
                    JOptionPane.showMessageDialog(SequenceDbDetailsDialog.this,
                            "File " + selectedFastaFile + " not found.",
                            "FASTA Import Error", JOptionPane.WARNING_MESSAGE);
                    e.printStackTrace();
                    return;
                } catch (Exception e) {
                    progressDialog.setRunFinished();
                    JOptionPane.showMessageDialog(SequenceDbDetailsDialog.this, JOptionEditorPane.getJOptionEditorPane(
                            "There was an error importing the FASTA file:<br>"
                            + e.getMessage() + "<br>"
                            + "See <a href=\"https://compomics.github.io/projects/searchgui/wiki/databasehelp.html\">DatabaseHelp</a> for help."),
                            "FASTA Import Error", JOptionPane.WARNING_MESSAGE);
                    e.printStackTrace();
                    return;
                }

                progressDialog.setRunFinished();

            }
        }.start();
    }

    /**
     * Appends decoy sequences to the given target database file.
     *
     * @param targetFile the target database file
     * @param progressDialog the progress dialog
     */
    private void generateTargetDecoyDatabase() {

        // set up the new fasta file name
        String newFasta = selectedFastaFile;
        File originalFastaFile = new File(selectedFastaFile);

        // remove the ending .fasta (if there)
        if (selectedFastaFile.lastIndexOf(".") != -1) {
            newFasta = selectedFastaFile.substring(0, selectedFastaFile.lastIndexOf("."));
        }

        // add the target decoy tag
        newFasta += fastaParameters.getTargetDecoyFileNameSuffix() + ".fasta";

        try {

            File newFile = new File(newFasta);

            progressDialog.setTitle("Appending Decoy Sequences. Please Wait...");
            progressDialog.setPrimaryProgressCounterIndeterminate(false);
            progressDialog.setPrimaryProgressCounter(0);
            progressDialog.setMaxPrimaryProgressCounter(fastaSummary.nSequences);

            DecoyConverter.appendDecoySequences(originalFastaFile, newFile, fastaParameters, progressDialog);

            progressDialog.setTitle("Getting Database Details. Please Wait...");

            progressDialog.setPrimaryProgressCounterIndeterminate(true);

            selectedFastaFile = newFile.getAbsolutePath();
            fastaParameters = DecoyConverter.getDecoyParameters(fastaParameters);
            fastaSummary = DecoyConverter.getDecoySummary(originalFastaFile, fastaSummary);

        } catch (OutOfMemoryError error) {

            Runtime.getRuntime().gc();
            JOptionPane.showMessageDialog(SequenceDbDetailsDialog.this,
                    "The tool used up all the available memory and had to be stopped.\n"
                    + "Memory boundaries are set in the Edit menu (Edit > Java Options).",
                    "Out Of Memory Error",
                    JOptionPane.ERROR_MESSAGE);
            System.out.println("Ran out of memory!");
            error.printStackTrace();

        } catch (FileNotFoundException e) {

            JOptionPane.showMessageDialog(SequenceDbDetailsDialog.this,
                    new String[]{"FASTA Import Error.", "File " + selectedFastaFile + " not found."},
                    "FASTA Import Error", JOptionPane.WARNING_MESSAGE);
            e.printStackTrace();

        } catch (Exception e) {

            JOptionPane.showMessageDialog(SequenceDbDetailsDialog.this,
                    new String[]{"FASTA Import Error.", "File " + selectedFastaFile + " could not be imported."},
                    "FASTA Import Error", JOptionPane.WARNING_MESSAGE);
            e.printStackTrace();

        }
    }

    /**
     * Copies the content of the FASTA file to a new file and replaces any white
     * space in the file name with '_' instead. Returns the new file, null if an
     * error occurred.
     *
     * @param file the FASTA file to rename
     * @return the renamed FASTA file
     */
    public File renameFastaFileName(File file) {

        String tempName = file.getName();
        tempName = tempName.replaceAll(" ", "_");

        File renamedFile = new File(file.getParentFile().getAbsolutePath() + File.separator + tempName);

        boolean success = false;

        try {

            success = renamedFile.createNewFile();

            if (success) {

                IoUtil.copyFile(file, renamedFile);

            }

        } catch (IOException e) {

            JOptionPane.showMessageDialog(this, "An error occurred while renaming the file.",
                    "Please Rename File", JOptionPane.WARNING_MESSAGE);
            e.printStackTrace();
            success = false;

        }

        if (success) {

            JOptionPane.showMessageDialog(this, "Your FASTA file name contained white space and has been renamed to:\n"
                    + file.getParentFile().getAbsolutePath() + File.separator + tempName, "Renamed File", JOptionPane.WARNING_MESSAGE);

            return renamedFile;

        }

        return null;
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        backgroundPanel = new javax.swing.JPanel();
        cancelButton = new javax.swing.JButton();
        okButton = new javax.swing.JButton();
        databaseInformationPanel = new javax.swing.JPanel();
        nameLabel = new javax.swing.JLabel();
        dbNameTxt = new javax.swing.JTextField();
        typeLabel = new javax.swing.JLabel();
        fileTxt = new javax.swing.JTextField();
        versionLabel = new javax.swing.JLabel();
        versionTxt = new javax.swing.JTextField();
        lastModifiedLabel = new javax.swing.JLabel();
        lastModifiedTxt = new javax.swing.JTextField();
        sizeLabel = new javax.swing.JLabel();
        sizeTxt = new javax.swing.JTextField();
        decoyButton = new javax.swing.JButton();
        browseButton = new javax.swing.JButton();
        fileLabel = new javax.swing.JLabel();
        typeJTextField = new javax.swing.JTextField();
        speciesJTextField = new javax.swing.JTextField();
        speciesLabel = new javax.swing.JLabel();
        previewPanel = new javax.swing.JPanel();
        proteinYxtScrollPane = new javax.swing.JScrollPane();
        proteinTxt = new javax.swing.JTextArea();
        proteinLabel = new javax.swing.JLabel();
        accessionsSpinner = new javax.swing.JSpinner();
        targetDecoyTxt = new javax.swing.JLabel();
        editFastaParametersJLabel = new javax.swing.JLabel();
        databaseHelpSettingsJLabel = new javax.swing.JLabel();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Database Overview");
        setMinimumSize(new java.awt.Dimension(500, 500));
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(java.awt.event.WindowEvent evt) {
                formWindowClosing(evt);
            }
        });

        backgroundPanel.setBackground(new java.awt.Color(230, 230, 230));

        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cancelButtonActionPerformed(evt);
            }
        });

        okButton.setText("OK");
        okButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                okButtonActionPerformed(evt);
            }
        });

        databaseInformationPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Database Details"));
        databaseInformationPanel.setOpaque(false);

        nameLabel.setText("Name");

        dbNameTxt.setEditable(false);
        dbNameTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);

        typeLabel.setText("Type(s)");

        fileTxt.setEditable(false);
        fileTxt.setHorizontalAlignment(javax.swing.JTextField.LEFT);

        versionLabel.setText("Version");

        versionTxt.setEditable(false);
        versionTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);

        lastModifiedLabel.setText("Modified");

        lastModifiedTxt.setEditable(false);
        lastModifiedTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);

        sizeLabel.setText("Size");

        sizeTxt.setEditable(false);
        sizeTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);

        decoyButton.setText("Decoys");
        decoyButton.setPreferredSize(new java.awt.Dimension(75, 25));
        decoyButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                decoyButtonActionPerformed(evt);
            }
        });

        browseButton.setText("Browse");
        browseButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                browseButtonActionPerformed(evt);
            }
        });

        fileLabel.setText("File");

        typeJTextField.setEditable(false);
        typeJTextField.setHorizontalAlignment(javax.swing.JTextField.CENTER);

        speciesJTextField.setEditable(false);
        speciesJTextField.setHorizontalAlignment(javax.swing.JTextField.CENTER);

        speciesLabel.setText("Species");

        javax.swing.GroupLayout databaseInformationPanelLayout = new javax.swing.GroupLayout(databaseInformationPanel);
        databaseInformationPanel.setLayout(databaseInformationPanelLayout);
        databaseInformationPanelLayout.setHorizontalGroup(
            databaseInformationPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(databaseInformationPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(databaseInformationPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, databaseInformationPanelLayout.createSequentialGroup()
                        .addComponent(fileLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 75, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fileTxt)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(browseButton, javax.swing.GroupLayout.PREFERRED_SIZE, 75, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(decoyButton, javax.swing.GroupLayout.PREFERRED_SIZE, 75, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(databaseInformationPanelLayout.createSequentialGroup()
                        .addComponent(sizeLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 75, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(sizeTxt))
                    .addGroup(databaseInformationPanelLayout.createSequentialGroup()
                        .addComponent(typeLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 75, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(typeJTextField))
                    .addGroup(databaseInformationPanelLayout.createSequentialGroup()
                        .addComponent(versionLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 75, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(versionTxt))
                    .addGroup(databaseInformationPanelLayout.createSequentialGroup()
                        .addComponent(nameLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 75, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(dbNameTxt))
                    .addGroup(databaseInformationPanelLayout.createSequentialGroup()
                        .addComponent(lastModifiedLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 75, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(lastModifiedTxt))
                    .addGroup(databaseInformationPanelLayout.createSequentialGroup()
                        .addComponent(speciesLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 75, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(speciesJTextField)))
                .addContainerGap())
        );

        databaseInformationPanelLayout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {fileLabel, lastModifiedLabel, nameLabel, sizeLabel, typeLabel, versionLabel});

        databaseInformationPanelLayout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {browseButton, decoyButton});

        databaseInformationPanelLayout.setVerticalGroup(
            databaseInformationPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(databaseInformationPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(databaseInformationPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fileTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(decoyButton, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(browseButton)
                    .addComponent(fileLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(databaseInformationPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(nameLabel)
                    .addComponent(dbNameTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(databaseInformationPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(speciesLabel)
                    .addComponent(speciesJTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(databaseInformationPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(typeLabel)
                    .addComponent(typeJTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(databaseInformationPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(versionTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(versionLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(databaseInformationPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(sizeLabel)
                    .addComponent(sizeTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(databaseInformationPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(lastModifiedLabel)
                    .addComponent(lastModifiedTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        databaseInformationPanelLayout.linkSize(javax.swing.SwingConstants.VERTICAL, new java.awt.Component[] {browseButton, decoyButton});

        previewPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Preview"));
        previewPanel.setOpaque(false);

        proteinTxt.setEditable(false);
        proteinTxt.setColumns(20);
        proteinTxt.setLineWrap(true);
        proteinTxt.setRows(5);
        proteinTxt.setWrapStyleWord(true);
        proteinYxtScrollPane.setViewportView(proteinTxt);

        proteinLabel.setText("Protein");

        accessionsSpinner.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                accessionsSpinnerStateChanged(evt);
            }
        });

        targetDecoyTxt.setText("(target/decoy)");

        editFastaParametersJLabel.setForeground(new java.awt.Color(0, 0, 255));
        editFastaParametersJLabel.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        editFastaParametersJLabel.setText("<html><a style=\"text-decoration: none\">Edit FASTA parsing options</a></html>");
        editFastaParametersJLabel.setToolTipText("Edit the FASTA parsing options");
        editFastaParametersJLabel.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                editFastaParametersJLabelMouseClicked(evt);
            }
            public void mouseEntered(java.awt.event.MouseEvent evt) {
                editFastaParametersJLabelMouseEntered(evt);
            }
            public void mouseExited(java.awt.event.MouseEvent evt) {
                editFastaParametersJLabelMouseExited(evt);
            }
        });

        javax.swing.GroupLayout previewPanelLayout = new javax.swing.GroupLayout(previewPanel);
        previewPanel.setLayout(previewPanelLayout);
        previewPanelLayout.setHorizontalGroup(
            previewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(previewPanelLayout.createSequentialGroup()
                .addGroup(previewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(previewPanelLayout.createSequentialGroup()
                        .addGap(16, 16, 16)
                        .addComponent(proteinLabel)
                        .addGap(18, 18, 18)
                        .addComponent(accessionsSpinner, javax.swing.GroupLayout.PREFERRED_SIZE, 192, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(targetDecoyTxt)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 199, Short.MAX_VALUE)
                        .addComponent(editFastaParametersJLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 155, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(proteinYxtScrollPane))
                .addContainerGap())
        );
        previewPanelLayout.setVerticalGroup(
            previewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(previewPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(proteinYxtScrollPane, javax.swing.GroupLayout.DEFAULT_SIZE, 160, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(previewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(proteinLabel)
                    .addComponent(accessionsSpinner, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(targetDecoyTxt)
                    .addComponent(editFastaParametersJLabel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap())
        );

        databaseHelpSettingsJLabel.setForeground(new java.awt.Color(0, 0, 255));
        databaseHelpSettingsJLabel.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        databaseHelpSettingsJLabel.setText("<html><a style=\"text-decoration: none\">Database help?</a></html>");
        databaseHelpSettingsJLabel.setToolTipText("Open Database Help");
        databaseHelpSettingsJLabel.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                databaseHelpSettingsJLabelMouseClicked(evt);
            }
            public void mouseEntered(java.awt.event.MouseEvent evt) {
                databaseHelpSettingsJLabelMouseEntered(evt);
            }
            public void mouseExited(java.awt.event.MouseEvent evt) {
                databaseHelpSettingsJLabelMouseExited(evt);
            }
        });

        javax.swing.GroupLayout backgroundPanelLayout = new javax.swing.GroupLayout(backgroundPanel);
        backgroundPanel.setLayout(backgroundPanelLayout);
        backgroundPanelLayout.setHorizontalGroup(
            backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(backgroundPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(backgroundPanelLayout.createSequentialGroup()
                        .addGap(10, 10, 10)
                        .addComponent(databaseHelpSettingsJLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(okButton, javax.swing.GroupLayout.PREFERRED_SIZE, 75, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(cancelButton, javax.swing.GroupLayout.PREFERRED_SIZE, 75, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(previewPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(databaseInformationPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap())
        );

        backgroundPanelLayout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {cancelButton, okButton});

        backgroundPanelLayout.setVerticalGroup(
            backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(backgroundPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(databaseInformationPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(previewPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(cancelButton)
                    .addComponent(okButton)
                    .addComponent(databaseHelpSettingsJLabel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap())
        );

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(backgroundPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(backgroundPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    /**
     * Saves changes and closes the dialog
     *
     * @param evt the action event
     */
    private void okButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_okButtonActionPerformed

        UtilitiesUserParameters.saveUserParameters(utilitiesUserParameters);

        dispose();

    }//GEN-LAST:event_okButtonActionPerformed

    /**
     * Close the dialog.
     *
     * @param evt the action event
     */
    private void cancelButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cancelButtonActionPerformed

        canceled = true;

        dispose();

    }//GEN-LAST:event_cancelButtonActionPerformed

    /**
     * Open a file chooser to select a FASTA file.
     *
     * @param evt the action event
     */
    private void browseButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_browseButtonActionPerformed
        selectDB(false);
    }//GEN-LAST:event_browseButtonActionPerformed

    /**
     * Add decoys.
     *
     * @param evt the action event
     */
    private void decoyButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_decoyButtonActionPerformed

        progressDialog = new ProgressDialogX(this, parentFrame,
                normalImange,
                waitingImage,
                true);
        progressDialog.setPrimaryProgressCounterIndeterminate(true);
        progressDialog.setTitle("Creating Decoy. Please Wait...");

        new Thread(new Runnable() {
            public void run() {
                try {
                    progressDialog.setVisible(true);
                } catch (IndexOutOfBoundsException e) {
                    // ignore
                }
            }
        }, "ProgressDialog").start();

        new Thread("DecoyThread") {
            public void run() {

                generateTargetDecoyDatabase();

                if (!progressDialog.isRunCanceled()) {

                    setUpGUI();

                }

                progressDialog.setRunFinished();

            }
        }.start();

    }//GEN-LAST:event_decoyButtonActionPerformed

    /**
     * Update the sequence.
     *
     * @param evt the change event
     */
    private void accessionsSpinnerStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_accessionsSpinnerStateChanged
        updateSequence();
    }//GEN-LAST:event_accessionsSpinnerStateChanged

    /**
     * Open the database help page.
     *
     * @param evt the mouse event
     */
    private void databaseHelpSettingsJLabelMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_databaseHelpSettingsJLabelMouseClicked
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        BareBonesBrowserLaunch.openURL("https://compomics.github.io/projects/searchgui/wiki/databasehelp.html");
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }//GEN-LAST:event_databaseHelpSettingsJLabelMouseClicked

    /**
     * Change the cursor to a hand cursor.
     *
     * @param evt the mouse event
     */
    private void databaseHelpSettingsJLabelMouseEntered(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_databaseHelpSettingsJLabelMouseEntered
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.HAND_CURSOR));
    }//GEN-LAST:event_databaseHelpSettingsJLabelMouseEntered

    /**
     * Change cursor back to the default cursor.
     *
     * @param evt the mouse event
     */
    private void databaseHelpSettingsJLabelMouseExited(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_databaseHelpSettingsJLabelMouseExited
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }//GEN-LAST:event_databaseHelpSettingsJLabelMouseExited

    /**
     * Close the dialog.
     *
     * @param evt the action event
     */
    private void formWindowClosing(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowClosing
        cancelButtonActionPerformed(null);
    }//GEN-LAST:event_formWindowClosing

    /**
     * Open the FASTA parameters dialog.
     *
     * @param evt
     */
    private void editFastaParametersJLabelMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_editFastaParametersJLabelMouseClicked
        FastaParametersDialog fastaParametersDialog = new FastaParametersDialog(this, parentFrame, fastaParameters, dbEditable);

        if (!fastaParametersDialog.isCanceled()) {
            FastaParameters newFastaParameters = fastaParametersDialog.getFastaSettings();

            if (!newFastaParameters.equals(fastaParameters)) {
                fastaParameters = newFastaParameters;
                loadFastaFile(selectedFastaFile, false, true, true);
            }
        }
    }//GEN-LAST:event_editFastaParametersJLabelMouseClicked

    /**
     * Change the cursor to a hand cursor.
     *
     * @param evt the mouse event
     */
    private void editFastaParametersJLabelMouseEntered(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_editFastaParametersJLabelMouseEntered
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.HAND_CURSOR));
    }//GEN-LAST:event_editFastaParametersJLabelMouseEntered

    /**
     * Change cursor back to the default cursor.
     *
     * @param evt the mouse event
     */
    private void editFastaParametersJLabelMouseExited(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_editFastaParametersJLabelMouseExited
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }//GEN-LAST:event_editFastaParametersJLabelMouseExited

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JSpinner accessionsSpinner;
    private javax.swing.JPanel backgroundPanel;
    private javax.swing.JButton browseButton;
    private javax.swing.JButton cancelButton;
    private javax.swing.JLabel databaseHelpSettingsJLabel;
    private javax.swing.JPanel databaseInformationPanel;
    private javax.swing.JTextField dbNameTxt;
    private javax.swing.JButton decoyButton;
    private javax.swing.JLabel editFastaParametersJLabel;
    private javax.swing.JLabel fileLabel;
    private javax.swing.JTextField fileTxt;
    private javax.swing.JLabel lastModifiedLabel;
    private javax.swing.JTextField lastModifiedTxt;
    private javax.swing.JLabel nameLabel;
    private javax.swing.JButton okButton;
    private javax.swing.JPanel previewPanel;
    private javax.swing.JLabel proteinLabel;
    private javax.swing.JTextArea proteinTxt;
    private javax.swing.JScrollPane proteinYxtScrollPane;
    private javax.swing.JLabel sizeLabel;
    private javax.swing.JTextField sizeTxt;
    private javax.swing.JTextField speciesJTextField;
    private javax.swing.JLabel speciesLabel;
    private javax.swing.JLabel targetDecoyTxt;
    private javax.swing.JTextField typeJTextField;
    private javax.swing.JLabel typeLabel;
    private javax.swing.JLabel versionLabel;
    private javax.swing.JTextField versionTxt;
    // End of variables declaration//GEN-END:variables

    /**
     * Buffers proteins sampled from the database.
     *
     * @throws IOException exception thrown if an error occurred while reading
     * the FASTA file
     */
    private void bufferProteins() throws IOException {

        int i = 0, previousSize = proteinsSample.size();

        Protein protein;

        while (i < SAMPLE_BATCH_SIZE && (protein = proteinIterator.getNextProtein()) != null) {

            String accession = protein.getAccession();
            accessionsSample.add(accession);
            proteinsSample.put(accession, protein);
            headersSample.put(accession, proteinIterator.getLastHeader());

        }

        accessionsSpinner.setModel(new SpinnerListModel(accessionsSample));
        accessionsSpinner.setValue(accessionsSample.get(previousSize));

    }

    /**
     * Returns the selected FASTA file.
     *
     * @return the selected FASTA file
     */
    public String getSelectedFastaFile() {
        return selectedFastaFile;
    }

    /**
     * Returns the FASTA parameters.
     *
     * @return the FASTA parameters
     */
    public FastaParameters getFastaParameters() {
        return fastaParameters;
    }

    /**
     * Returns a boolean indicating whether the database selection was canceled
     * by the user.
     *
     * @return a boolean indicating whether the database selection was canceled
     * by the user
     */
    public boolean isCanceled() {
        return canceled;
    }
}
