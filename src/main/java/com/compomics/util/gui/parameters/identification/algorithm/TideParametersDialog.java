package com.compomics.util.gui.parameters.identification.algorithm;

import com.compomics.util.examples.BareBonesBrowserLaunch;
import com.compomics.util.gui.parameters.identification.IdentificationAlgorithmParameter;
import com.compomics.util.parameters.identification.tool_specific.TideParameters;
import com.compomics.util.gui.GuiUtilities;
import com.compomics.util.gui.JOptionEditorPane;
import java.awt.Color;
import java.awt.Dialog;
import javax.swing.JOptionPane;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import com.compomics.util.gui.parameters.identification.AlgorithmParametersDialog;

/**
 * Dialog for the Tide specific settings.
 *
 * @author Harald Barsnes
 */
public class TideParametersDialog extends javax.swing.JDialog implements AlgorithmParametersDialog {

    /**
     * Boolean indicating whether the used canceled the editing.
     */
    private boolean cancelled = false;
    /**
     * Boolean indicating whether the settings can be edited by the user.
     */
    private boolean editable;

    /**
     * Creates a new TideSettingsDialog with a frame as owner.
     *
     * @param parent the parent frame
     * @param tideParameters the Tide parameters
     * @param editable boolean indicating whether the settings can be edited by
     * the user
     */
    public TideParametersDialog(java.awt.Frame parent, TideParameters tideParameters, boolean editable) {
        super(parent, true);
        this.editable = editable;
        initComponents();
        setUpGUI();
        populateGUI(tideParameters);
        validateInput(false);
        setLocationRelativeTo(parent);
        setVisible(true);
    }

    /**
     * Creates a new TideSettingsDialog with a dialog as owner.
     *
     * @param owner the dialog owner
     * @param parent the parent frame
     * @param tideParameters the Tide parameters
     * @param editable boolean indicating whether the settings can be edited by
     * the user
     */
    public TideParametersDialog(Dialog owner, java.awt.Frame parent, TideParameters tideParameters, boolean editable) {
        super(owner, true);
        this.editable = editable;
        initComponents();
        setUpGUI();
        populateGUI(tideParameters);
        validateInput(false);
        setLocationRelativeTo(owner);
        setVisible(true);
    }

    /**
     * Sets up the GUI.
     */
    private void setUpGUI() {

        removePrecursorPeakCombo.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        enzymeTypeCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        useFlankingCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        removePrecursorPeakCombo.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        monoPrecursorCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        peptideListCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        decoyFormatCombo.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        keepTerminalAaCombo.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        removeMethionineCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        exactPvalueCombo.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        scoreFunctionCombo.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        useNeutralLossCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        outputFormatCombo.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        removeTempFoldersCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));

        minPepLengthTxt.setEditable(editable);
        minPepLengthTxt.setEnabled(editable);
        maxPepLengthTxt.setEditable(editable);
        maxPepLengthTxt.setEnabled(editable);
        minPrecursorMassTxt.setEditable(editable);
        minPrecursorMassTxt.setEnabled(editable);
        maxPrecursorMassTxt.setEditable(editable);
        maxPrecursorMassTxt.setEnabled(editable);
        monoPrecursorCmb.setEnabled(editable);
        removeMethionineCmb.setEnabled(editable);
        //minPtmsPerPeptideTxt.setEditable(editable);
        //minPtmsPerPeptideTxt.setEnabled(editable);
        maxPtmsPerPeptideTxt.setEditable(editable);
        maxPtmsPerPeptideTxt.setEnabled(editable);
        maxVariablePtmsPerTypeTxt.setEditable(editable);
        maxVariablePtmsPerTypeTxt.setEnabled(editable);
        enzymeTypeCmb.setEnabled(editable);
        peptideListCmb.setEnabled(editable);
        decoyFormatCombo.setEnabled(editable);
        keepTerminalAaCombo.setEnabled(editable);
        decoySeedTxt.setEditable(editable);
        decoySeedTxt.setEnabled(editable);
        removeTempFoldersCmb.setEnabled(editable);
        exactPvalueCombo.setEnabled(editable);
        scoreFunctionCombo.setEnabled(editable);
        minSpectrumMzTxt.setEditable(editable);
        minSpectrumMzTxt.setEnabled(editable);
        maxSpectrumMzTxt.setEditable(editable);
        maxSpectrumMzTxt.setEnabled(editable);
        minPeaksTxt.setEditable(editable);
        minPeaksTxt.setEnabled(editable);
        removePrecursorPeakCombo.setEnabled(editable);
        removePrecursorPeakToleranceTxt.setEditable(editable);
        removePrecursorPeakToleranceTxt.setEnabled(editable);
        useFlankingCmb.setEnabled(editable);
        useNeutralLossCmb.setEnabled(editable);
        mzBinWidthTxt.setEditable(editable);
        mzBinWidthTxt.setEnabled(editable);
        mzBinOffsetTxt.setEditable(editable);
        mzBinOffsetTxt.setEnabled(editable);
        numberMatchesTxt.setEditable(editable);
        numberMatchesTxt.setEnabled(editable);
        outputFormatCombo.setEnabled(editable);

    }

    /**
     * Populates the GUI using the given settings.
     *
     * @param tideParameters the parameters to display
     */
    private void populateGUI(TideParameters tideParameters) {

        if (tideParameters.getMinPeptideLength() != null) {
            minPepLengthTxt.setText(tideParameters.getMinPeptideLength() + "");
        }
        if (tideParameters.getMaxPeptideLength() != null) {
            maxPepLengthTxt.setText(tideParameters.getMaxPeptideLength() + "");
        }
        if (tideParameters.getMinPrecursorMass() != null) {
            minPrecursorMassTxt.setText(tideParameters.getMinPrecursorMass() + "");
        }
        if (tideParameters.getMaxPrecursorMass() != null) {
            maxPrecursorMassTxt.setText(tideParameters.getMaxPrecursorMass() + "");
        }
        if (tideParameters.getMonoisotopicPrecursor() != null) {
            if (tideParameters.getMonoisotopicPrecursor()) {
                monoPrecursorCmb.setSelectedIndex(0);
            } else {
                monoPrecursorCmb.setSelectedIndex(1);
            }
        }
        if (tideParameters.getClipNtermMethionine() != null) {
            if (tideParameters.getClipNtermMethionine()) {
                removeMethionineCmb.setSelectedIndex(0);
            } else {
                removeMethionineCmb.setSelectedIndex(1);
            }
        }
//        if (tideParameters.getMinVariableModificationsPerPeptide() != null) {
//            minPtmsPerPeptideTxt.setText(tideParameters.getMinVariableModificationsPerPeptide() + "");
//        }
        if (tideParameters.getMaxVariableModificationsPerPeptide() != null) {
            maxPtmsPerPeptideTxt.setText(tideParameters.getMaxVariableModificationsPerPeptide() + "");
        }
        if (tideParameters.getMaxVariableModificationsPerTypePerPeptide() != null) {
            maxVariablePtmsPerTypeTxt.setText(tideParameters.getMaxVariableModificationsPerTypePerPeptide() + "");
        }
        if (tideParameters.getDigestionType() != null) {
            enzymeTypeCmb.setSelectedItem(tideParameters.getDigestionType());
        }
        if (tideParameters.getPrintPeptides() != null) {
            if (tideParameters.getPrintPeptides()) {
                peptideListCmb.setSelectedIndex(0);
            } else {
                peptideListCmb.setSelectedIndex(1);
            }
        }
        if (tideParameters.getDecoyFormat() != null) {
            decoyFormatCombo.setSelectedItem(tideParameters.getDecoyFormat());
        }
        if (tideParameters.getKeepTerminalAminoAcids() != null) {
            keepTerminalAaCombo.setSelectedItem(tideParameters.getKeepTerminalAminoAcids());
        }
        if (tideParameters.getDecoySeed() != null) {
            decoySeedTxt.setText(tideParameters.getDecoySeed() + "");
        }
        if (tideParameters.getRemoveTempFolders() != null) {
            if (tideParameters.getRemoveTempFolders()) {
                removeTempFoldersCmb.setSelectedIndex(0);
            } else {
                removeTempFoldersCmb.setSelectedIndex(1);
            }
        }
        if (tideParameters.getComputeExactPValues() != null) {
            if (tideParameters.getComputeExactPValues()) {
                exactPvalueCombo.setSelectedIndex(0);
            } else {
                exactPvalueCombo.setSelectedIndex(1);
            }
        }
        if (tideParameters.getScoreFunction() != null) {
            scoreFunctionCombo.setSelectedItem(tideParameters.getScoreFunction());
        }
        if (tideParameters.getMinSpectrumMz() != null) {
            minSpectrumMzTxt.setText(tideParameters.getMinSpectrumMz() + "");
        }
        if (tideParameters.getMaxSpectrumMz() != null) {
            maxSpectrumMzTxt.setText(tideParameters.getMaxSpectrumMz() + "");
        }
        if (tideParameters.getMinSpectrumPeaks() != null) {
            minPeaksTxt.setText(tideParameters.getMinSpectrumPeaks() + "");
        }
        if (tideParameters.getRemovePrecursor() != null) {
            if (tideParameters.getRemovePrecursor()) {
                removePrecursorPeakCombo.setSelectedIndex(0);
            } else {
                removePrecursorPeakCombo.setSelectedIndex(1);
            }
        }
        if (tideParameters.getRemovePrecursorTolerance() != null) {
            removePrecursorPeakToleranceTxt.setText("" + tideParameters.getRemovePrecursorTolerance());
        }
        if (tideParameters.getUseFlankingPeaks() != null) {
            if (tideParameters.getUseFlankingPeaks()) {
                useFlankingCmb.setSelectedIndex(0);
            } else {
                useFlankingCmb.setSelectedIndex(1);
            }
        }
        if (tideParameters.getUseNeutralLossPeaks() != null) {
            if (tideParameters.getUseNeutralLossPeaks()) {
                useNeutralLossCmb.setSelectedIndex(0);
            } else {
                useNeutralLossCmb.setSelectedIndex(1);
            }
        }
        if (tideParameters.getMzBinWidth() != null) {
            mzBinWidthTxt.setText("" + tideParameters.getMzBinWidth());
        }
        if (tideParameters.getMzBinOffset() != null) {
            mzBinOffsetTxt.setText("" + tideParameters.getMzBinOffset());
        }
        if (tideParameters.getNumberOfSpectrumMatches() != null) {
            numberMatchesTxt.setText("" + tideParameters.getNumberOfSpectrumMatches());
        }
        if (tideParameters.getTextOutput()) {
            outputFormatCombo.setSelectedItem("Text");
        } else if (tideParameters.getSqtOutput()) {
            outputFormatCombo.setSelectedItem("SQT");
        } else if (tideParameters.getPepXmlOutput()) {
            outputFormatCombo.setSelectedItem("pepxml");
        } else if (tideParameters.getMzidOutput()) {
            outputFormatCombo.setSelectedItem("mzIdentML");
        } else if (tideParameters.getPinOutput()) {
            outputFormatCombo.setSelectedItem("Percolator input file");
        }
    }

    @Override
    public boolean isCancelled() {
        return cancelled;
    }

    @Override
    public IdentificationAlgorithmParameter getParameters() {
        return getInput();
    }

    /**
     * Returns the user selection as Tide parameters object.
     *
     * @return the user selection
     */
    public TideParameters getInput() {

        TideParameters result = new TideParameters();

        String input = minPepLengthTxt.getText().trim();
        if (!input.equals("")) {
            result.setMinPeptideLength(Integer.valueOf(input));
        }
        input = maxPepLengthTxt.getText().trim();
        if (!input.equals("")) {
            result.setMaxPeptideLength(Integer.valueOf(input));
        }

        input = minPrecursorMassTxt.getText().trim();
        if (!input.equals("")) {
            result.setMinPrecursorMass(Double.valueOf(input));
        }
        input = maxPrecursorMassTxt.getText().trim();
        if (!input.equals("")) {
            result.setMaxPrecursorMass(Double.valueOf(input));
        }

        result.setMonoisotopicPrecursor(monoPrecursorCmb.getSelectedIndex() == 0);
        result.setClipNtermMethionine(removeMethionineCmb.getSelectedIndex() == 0);

//        input = minPtmsPerPeptideTxt.getText().trim();
//        if (!input.equals("")) {
//            result.setMinVariableModificationsPerPeptide(Integer.valueOf(input));
//        }
        input = maxPtmsPerPeptideTxt.getText().trim();
        if (!input.equals("")) {
            result.setMaxVariableModificationsPerPeptide(Integer.valueOf(input));
        }

        input = maxVariablePtmsPerTypeTxt.getText().trim();
        if (!input.equals("")) {
            result.setMaxVariableModificationsPerTypePerPeptide(Integer.valueOf(input));
        }

        result.setDigestionType((String) enzymeTypeCmb.getSelectedItem());
        result.setPrintPeptides(peptideListCmb.getSelectedIndex() == 0);
        result.setDecoyFormat((String) decoyFormatCombo.getSelectedItem());
        result.setRemoveTempFolders(removeTempFoldersCmb.getSelectedIndex() == 0);
        result.setKeepTerminalAminoAcids((String) keepTerminalAaCombo.getSelectedItem());

        input = decoySeedTxt.getText().trim();
        if (!input.equals("")) {
            result.setDecoySeed(Integer.valueOf(input));
        }

        result.setComputeExactPValues(exactPvalueCombo.getSelectedIndex() == 0);
        result.setScoreFunction((String) scoreFunctionCombo.getSelectedItem());

        input = minSpectrumMzTxt.getText().trim();
        if (!input.equals("")) {
            result.setMinSpectrumMz(Double.valueOf(input));
        }
        input = maxSpectrumMzTxt.getText().trim();
        if (!input.equals("")) {
            result.setMaxSpectrumMz(Double.valueOf(input));
        }

        input = minPeaksTxt.getText().trim();
        if (!input.equals("")) {
            result.setMinSpectrumPeaks(Integer.valueOf(input));
        }

        result.setRemovePrecursor(removePrecursorPeakCombo.getSelectedIndex() == 0);

        input = removePrecursorPeakToleranceTxt.getText().trim();
        if (!input.equals("")) {
            result.setRemovePrecursorTolerance(Double.valueOf(input));
        }

        result.setUseFlankingPeaks(useFlankingCmb.getSelectedIndex() == 0);
        result.setUseNeutralLossPeaks(useNeutralLossCmb.getSelectedIndex() == 0);

        input = mzBinWidthTxt.getText().trim();
        if (!input.equals("")) {
            result.setMzBinWidth(Double.valueOf(input));
        }

        input = mzBinOffsetTxt.getText().trim();
        if (!input.equals("")) {
            result.setMzBinOffset(Double.valueOf(input));
        }

        input = numberMatchesTxt.getText().trim();
        if (!input.equals("")) {
            result.setNumberOfSpectrumMatches(Integer.valueOf(input));
        }

        int selectedIndex = outputFormatCombo.getSelectedIndex();
        result.setTextOutput(selectedIndex == 0);
        result.setSqtOutput(selectedIndex == 1);
        result.setPepXmlOutput(selectedIndex == 2);
        result.setMzidOutput(selectedIndex == 3);
        result.setPinOutput(selectedIndex == 4);

        return result;
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
        tabbedPane = new javax.swing.JTabbedPane();
        indexPanel = new javax.swing.JPanel();
        peptideLengthLabel = new javax.swing.JLabel();
        minPepLengthTxt = new javax.swing.JTextField();
        peptideLengthDividerLabel = new javax.swing.JLabel();
        maxPepLengthTxt = new javax.swing.JTextField();
        precursorMassLabel = new javax.swing.JLabel();
        minPrecursorMassTxt = new javax.swing.JTextField();
        precursorMassDividerLabel = new javax.swing.JLabel();
        maxPrecursorMassTxt = new javax.swing.JTextField();
        monoPrecursorLabel = new javax.swing.JLabel();
        monoPrecursorCmb = new javax.swing.JComboBox();
        removeMethionineLabel = new javax.swing.JLabel();
        removeMethionineCmb = new javax.swing.JComboBox();
        variablePtmsPerPeptideLabel = new javax.swing.JLabel();
        maxPtmsPerPeptideTxt = new javax.swing.JTextField();
        decoyformatLabel = new javax.swing.JLabel();
        decoyFormatCombo = new javax.swing.JComboBox();
        keepTerminalAaLabel = new javax.swing.JLabel();
        keepTerminalAaCombo = new javax.swing.JComboBox();
        decoySeedLabel = new javax.swing.JLabel();
        decoySeedTxt = new javax.swing.JTextField();
        enzymeTypeLabel = new javax.swing.JLabel();
        enzymeTypeCmb = new javax.swing.JComboBox();
        peptideListLabel = new javax.swing.JLabel();
        peptideListCmb = new javax.swing.JComboBox();
        maxVariablePtmsPerTypeLabel = new javax.swing.JLabel();
        maxVariablePtmsPerTypeTxt = new javax.swing.JTextField();
        removeTempFoldersLabel = new javax.swing.JLabel();
        removeTempFoldersCmb = new javax.swing.JComboBox();
        minPtmsPerPeptideTxt = new javax.swing.JTextField();
        variablePtmsPerPeptideDividerLabel = new javax.swing.JLabel();
        searchPanel = new javax.swing.JPanel();
        useFlankingLabel = new javax.swing.JLabel();
        useFlankingCmb = new javax.swing.JComboBox();
        spectrumMzLabel = new javax.swing.JLabel();
        minSpectrumMzTxt = new javax.swing.JTextField();
        spectrumMzDividerLabel = new javax.swing.JLabel();
        maxSpectrumMzTxt = new javax.swing.JTextField();
        numberMatchesLabel = new javax.swing.JLabel();
        numberMatchesTxt = new javax.swing.JTextField();
        removePrecursorPeakLabel = new javax.swing.JLabel();
        removePrecursorPeakCombo = new javax.swing.JComboBox();
        removePrecursorPeakToleranceLbl = new javax.swing.JLabel();
        removePrecursorPeakToleranceTxt = new javax.swing.JTextField();
        mzBinWidthLabel = new javax.swing.JLabel();
        mzBinWidthTxt = new javax.swing.JTextField();
        minPeaksLbl = new javax.swing.JLabel();
        minPeaksTxt = new javax.swing.JTextField();
        exactPvalueLabel = new javax.swing.JLabel();
        exactPvalueCombo = new javax.swing.JComboBox();
        scoreFunctionLabel = new javax.swing.JLabel();
        scoreFunctionCombo = new javax.swing.JComboBox();
        useNeutralLossLabel = new javax.swing.JLabel();
        useNeutralLossCmb = new javax.swing.JComboBox();
        mzBinOffsetLabel = new javax.swing.JLabel();
        mzBinOffsetTxt = new javax.swing.JTextField();
        outputFormatLabel = new javax.swing.JLabel();
        outputFormatCombo = new javax.swing.JComboBox();
        openDialogHelpJButton = new javax.swing.JButton();
        advancedSettingsWarningLabel = new javax.swing.JLabel();
        okButton = new javax.swing.JButton();
        closeButton = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Tide Advanced Settings");
        setResizable(false);

        backgroundPanel.setBackground(new java.awt.Color(230, 230, 230));

        tabbedPane.setBackground(new java.awt.Color(230, 230, 230));
        tabbedPane.setOpaque(true);

        indexPanel.setBackground(new java.awt.Color(230, 230, 230));
        indexPanel.setPreferredSize(new java.awt.Dimension(518, 143));

        peptideLengthLabel.setText("Peptide Length (min - max)");

        minPepLengthTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        minPepLengthTxt.setText("6");
        minPepLengthTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                minPepLengthTxtKeyReleased(evt);
            }
        });

        peptideLengthDividerLabel.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        peptideLengthDividerLabel.setText("-");

        maxPepLengthTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        maxPepLengthTxt.setText("30");
        maxPepLengthTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                maxPepLengthTxtKeyReleased(evt);
            }
        });

        precursorMassLabel.setText("Precursor Mass (min - max)");

        minPrecursorMassTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        minPrecursorMassTxt.setText("0");
        minPrecursorMassTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                minPrecursorMassTxtKeyReleased(evt);
            }
        });

        precursorMassDividerLabel.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        precursorMassDividerLabel.setText("-");

        maxPrecursorMassTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        maxPrecursorMassTxt.setText("10000");
        maxPrecursorMassTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                maxPrecursorMassTxtKeyReleased(evt);
            }
        });

        monoPrecursorLabel.setText("Monoisotopic Precursor");

        monoPrecursorCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));

        removeMethionineLabel.setText("Remove Starting Methionine");

        removeMethionineCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));
        removeMethionineCmb.setSelectedIndex(1);

        variablePtmsPerPeptideLabel.setText("Variable PTMs per Peptide");

        maxPtmsPerPeptideTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        maxPtmsPerPeptideTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                maxPtmsPerPeptideTxtKeyReleased(evt);
            }
        });

        decoyformatLabel.setText("Decoy Format");

        decoyFormatCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "none", "shuffle", "peptide-reverse", "protein-reverse" }));

        keepTerminalAaLabel.setText("Keep Terminal AAs");

        keepTerminalAaCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "N", "C", "NC", "none" }));
        keepTerminalAaCombo.setSelectedIndex(2);

        decoySeedLabel.setText("Decoy Seed");

        decoySeedTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        decoySeedTxt.setText("1");
        decoySeedTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                decoySeedTxtKeyReleased(evt);
            }
        });

        enzymeTypeLabel.setText("Enzyme Type");

        enzymeTypeCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "full-digest", "partial-digest" }));

        peptideListLabel.setText("Peptide List");

        peptideListCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));
        peptideListCmb.setSelectedIndex(1);

        maxVariablePtmsPerTypeLabel.setText("Max Variable PTMs Per Type");

        maxVariablePtmsPerTypeTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        maxVariablePtmsPerTypeTxt.setText("2");
        maxVariablePtmsPerTypeTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                maxVariablePtmsPerTypeTxtKeyReleased(evt);
            }
        });

        removeTempFoldersLabel.setText("Recreate FASTA Index");

        removeTempFoldersCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));

        minPtmsPerPeptideTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        minPtmsPerPeptideTxt.setText("0");
        minPtmsPerPeptideTxt.setEnabled(false);
        minPtmsPerPeptideTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                minPtmsPerPeptideTxtKeyReleased(evt);
            }
        });

        variablePtmsPerPeptideDividerLabel.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        variablePtmsPerPeptideDividerLabel.setText("-");

        javax.swing.GroupLayout indexPanelLayout = new javax.swing.GroupLayout(indexPanel);
        indexPanel.setLayout(indexPanelLayout);
        indexPanelLayout.setHorizontalGroup(
            indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(indexPanelLayout.createSequentialGroup()
                .addGap(25, 25, 25)
                .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(indexPanelLayout.createSequentialGroup()
                        .addComponent(maxVariablePtmsPerTypeLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(maxVariablePtmsPerTypeTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(indexPanelLayout.createSequentialGroup()
                        .addComponent(peptideLengthLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(minPepLengthTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 71, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(12, 12, 12)
                        .addComponent(peptideLengthDividerLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 27, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(maxPepLengthTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 77, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(indexPanelLayout.createSequentialGroup()
                        .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addGroup(indexPanelLayout.createSequentialGroup()
                                .addComponent(precursorMassLabel, javax.swing.GroupLayout.DEFAULT_SIZE, 200, Short.MAX_VALUE)
                                .addGap(18, 18, 18))
                            .addGroup(indexPanelLayout.createSequentialGroup()
                                .addComponent(monoPrecursorLabel, javax.swing.GroupLayout.DEFAULT_SIZE, 200, Short.MAX_VALUE)
                                .addGap(18, 18, 18)))
                        .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(monoPrecursorCmb, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addGroup(indexPanelLayout.createSequentialGroup()
                                .addComponent(minPrecursorMassTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 71, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addGap(12, 12, 12)
                                .addComponent(precursorMassDividerLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 27, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addComponent(maxPrecursorMassTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 77, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addComponent(removeMethionineCmb, 0, 200, Short.MAX_VALUE)))
                    .addGroup(indexPanelLayout.createSequentialGroup()
                        .addComponent(keepTerminalAaLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(keepTerminalAaCombo, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(indexPanelLayout.createSequentialGroup()
                        .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(enzymeTypeLabel, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, 200, Short.MAX_VALUE)
                            .addComponent(peptideListLabel, javax.swing.GroupLayout.DEFAULT_SIZE, 200, Short.MAX_VALUE))
                        .addGap(18, 18, 18)
                        .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(peptideListCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(enzymeTypeCmb, 0, 200, Short.MAX_VALUE)))
                    .addGroup(indexPanelLayout.createSequentialGroup()
                        .addComponent(decoyformatLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(decoyFormatCombo, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(indexPanelLayout.createSequentialGroup()
                        .addComponent(variablePtmsPerPeptideLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(19, 19, 19)
                        .addComponent(minPtmsPerPeptideTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 77, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(variablePtmsPerPeptideDividerLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 39, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(maxPtmsPerPeptideTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 71, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(removeMethionineLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(indexPanelLayout.createSequentialGroup()
                        .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                            .addComponent(removeTempFoldersLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(decoySeedLabel, javax.swing.GroupLayout.DEFAULT_SIZE, 200, Short.MAX_VALUE))
                        .addGap(18, 18, 18)
                        .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(decoySeedTxt, javax.swing.GroupLayout.DEFAULT_SIZE, 200, Short.MAX_VALUE)
                            .addComponent(removeTempFoldersCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))))
                .addGap(0, 28, Short.MAX_VALUE))
        );

        indexPanelLayout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {enzymeTypeCmb, removeMethionineCmb});

        indexPanelLayout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {maxPrecursorMassTxt, minPrecursorMassTxt});

        indexPanelLayout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {maxPepLengthTxt, minPepLengthTxt});

        indexPanelLayout.setVerticalGroup(
            indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(indexPanelLayout.createSequentialGroup()
                .addGap(25, 25, 25)
                .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(minPepLengthTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(maxPepLengthTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(peptideLengthDividerLabel)
                    .addComponent(peptideLengthLabel))
                .addGap(0, 0, 0)
                .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(minPrecursorMassTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(maxPrecursorMassTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(precursorMassDividerLabel)
                    .addComponent(precursorMassLabel))
                .addGap(0, 0, 0)
                .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(monoPrecursorLabel)
                    .addComponent(monoPrecursorCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(removeMethionineLabel)
                    .addComponent(removeMethionineCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(minPtmsPerPeptideTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(variablePtmsPerPeptideDividerLabel))
                    .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(variablePtmsPerPeptideLabel)
                        .addComponent(maxPtmsPerPeptideTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addGap(0, 0, 0)
                .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(maxVariablePtmsPerTypeLabel)
                    .addComponent(maxVariablePtmsPerTypeTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(enzymeTypeLabel)
                    .addComponent(enzymeTypeCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(peptideListLabel)
                    .addComponent(peptideListCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(decoyformatLabel)
                    .addComponent(decoyFormatCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(keepTerminalAaLabel)
                    .addComponent(keepTerminalAaCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(decoySeedLabel)
                    .addComponent(decoySeedTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(indexPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(removeTempFoldersLabel)
                    .addComponent(removeTempFoldersCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap())
        );

        tabbedPane.addTab("Index", indexPanel);

        searchPanel.setBackground(new java.awt.Color(230, 230, 230));
        searchPanel.setPreferredSize(new java.awt.Dimension(518, 143));

        useFlankingLabel.setText("Use Flanking Peaks");

        useFlankingCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));
        useFlankingCmb.setSelectedIndex(1);
        useFlankingCmb.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                useFlankingCmbActionPerformed(evt);
            }
        });

        spectrumMzLabel.setText("Spectrum m/z (min - max)");

        minSpectrumMzTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        minSpectrumMzTxt.setText("0");
        minSpectrumMzTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                minSpectrumMzTxtKeyReleased(evt);
            }
        });

        spectrumMzDividerLabel.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        spectrumMzDividerLabel.setText("-");

        maxSpectrumMzTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        maxSpectrumMzTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                maxSpectrumMzTxtKeyReleased(evt);
            }
        });

        numberMatchesLabel.setText("Number of Spectrum Matches");

        numberMatchesTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        numberMatchesTxt.setText("10");
        numberMatchesTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                numberMatchesTxtKeyReleased(evt);
            }
        });

        removePrecursorPeakLabel.setText("Remove Precursor Peak (PP)");

        removePrecursorPeakCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));

        removePrecursorPeakToleranceLbl.setText("Remove PP Tolerance (in Da)");

        removePrecursorPeakToleranceTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        removePrecursorPeakToleranceTxt.setText("0.0");
        removePrecursorPeakToleranceTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                removePrecursorPeakToleranceTxtKeyReleased(evt);
            }
        });

        mzBinWidthLabel.setText("m/z Bin Width");

        mzBinWidthTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        mzBinWidthTxt.setText("0.02");
        mzBinWidthTxt.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                mzBinWidthTxtActionPerformed(evt);
            }
        });
        mzBinWidthTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                mzBinWidthTxtKeyReleased(evt);
            }
        });

        minPeaksLbl.setText("Minimum Number of Peaks");

        minPeaksTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        minPeaksTxt.setText("10");
        minPeaksTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                minPeaksTxtKeyReleased(evt);
            }
        });

        exactPvalueLabel.setText("Calculate Exact p-value");

        exactPvalueCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));

        scoreFunctionLabel.setText("Score Function");

        scoreFunctionCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "xcorr", "combined-p-values" }));

        useNeutralLossLabel.setText("Use Neutral Loss Peaks");

        useNeutralLossCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));
        useNeutralLossCmb.setSelectedIndex(1);

        mzBinOffsetLabel.setText("m/z Bin Offset");

        mzBinOffsetTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        mzBinOffsetTxt.setText("0.0");
        mzBinOffsetTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                mzBinOffsetTxtKeyReleased(evt);
            }
        });

        outputFormatLabel.setText("Output Format");

        outputFormatCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Text", "SQT", "pepxml", "mzIdentML", "Percolator input file" }));
        outputFormatCombo.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                outputFormatComboActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout searchPanelLayout = new javax.swing.GroupLayout(searchPanel);
        searchPanel.setLayout(searchPanelLayout);
        searchPanelLayout.setHorizontalGroup(
            searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(searchPanelLayout.createSequentialGroup()
                .addGap(25, 25, 25)
                .addGroup(searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(searchPanelLayout.createSequentialGroup()
                        .addComponent(outputFormatLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(outputFormatCombo, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(searchPanelLayout.createSequentialGroup()
                        .addComponent(mzBinOffsetLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(mzBinOffsetTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(searchPanelLayout.createSequentialGroup()
                        .addComponent(mzBinWidthLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(mzBinWidthTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(searchPanelLayout.createSequentialGroup()
                        .addComponent(minPeaksLbl, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(minPeaksTxt, javax.swing.GroupLayout.DEFAULT_SIZE, 200, Short.MAX_VALUE))
                    .addGroup(searchPanelLayout.createSequentialGroup()
                        .addGroup(searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                            .addComponent(spectrumMzLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(scoreFunctionLabel, javax.swing.GroupLayout.DEFAULT_SIZE, 200, Short.MAX_VALUE))
                        .addGap(18, 18, 18)
                        .addGroup(searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(scoreFunctionCombo, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addGroup(searchPanelLayout.createSequentialGroup()
                                .addComponent(minSpectrumMzTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 77, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(spectrumMzDividerLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 27, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addGap(18, 18, 18)
                                .addComponent(maxSpectrumMzTxt))))
                    .addGroup(searchPanelLayout.createSequentialGroup()
                        .addComponent(exactPvalueLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(exactPvalueCombo, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(searchPanelLayout.createSequentialGroup()
                        .addComponent(numberMatchesLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(numberMatchesTxt))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, searchPanelLayout.createSequentialGroup()
                        .addComponent(useNeutralLossLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(useNeutralLossCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addGroup(searchPanelLayout.createSequentialGroup()
                        .addGroup(searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(removePrecursorPeakLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(removePrecursorPeakToleranceLbl, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGap(18, 18, 18)
                        .addGroup(searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(removePrecursorPeakCombo, 0, 200, Short.MAX_VALUE)
                            .addComponent(removePrecursorPeakToleranceTxt)))
                    .addGroup(searchPanelLayout.createSequentialGroup()
                        .addComponent(useFlankingLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(useFlankingCmb, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(33, Short.MAX_VALUE))
        );

        searchPanelLayout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {minPeaksTxt, removePrecursorPeakCombo});

        searchPanelLayout.setVerticalGroup(
            searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(searchPanelLayout.createSequentialGroup()
                .addGap(25, 25, 25)
                .addGroup(searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(exactPvalueLabel)
                    .addComponent(exactPvalueCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(scoreFunctionLabel)
                    .addComponent(scoreFunctionCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(minSpectrumMzTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(maxSpectrumMzTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(spectrumMzDividerLabel)
                    .addComponent(spectrumMzLabel))
                .addGap(0, 0, 0)
                .addGroup(searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(minPeaksLbl)
                    .addComponent(minPeaksTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(removePrecursorPeakLabel)
                    .addComponent(removePrecursorPeakCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(removePrecursorPeakToleranceLbl)
                    .addComponent(removePrecursorPeakToleranceTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(useFlankingLabel)
                    .addComponent(useFlankingCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(useNeutralLossLabel)
                    .addComponent(useNeutralLossCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(mzBinWidthLabel)
                    .addComponent(mzBinWidthTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(mzBinOffsetLabel)
                    .addComponent(mzBinOffsetTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(numberMatchesLabel)
                    .addComponent(numberMatchesTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(searchPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(outputFormatLabel)
                    .addComponent(outputFormatCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(47, Short.MAX_VALUE))
        );

        tabbedPane.addTab("Search", searchPanel);

        openDialogHelpJButton.setIcon(new javax.swing.ImageIcon(getClass().getResource("/icons/help.GIF"))); // NOI18N
        openDialogHelpJButton.setToolTipText("Help");
        openDialogHelpJButton.setBorder(null);
        openDialogHelpJButton.setBorderPainted(false);
        openDialogHelpJButton.setContentAreaFilled(false);
        openDialogHelpJButton.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseEntered(java.awt.event.MouseEvent evt) {
                openDialogHelpJButtonMouseEntered(evt);
            }
            public void mouseExited(java.awt.event.MouseEvent evt) {
                openDialogHelpJButtonMouseExited(evt);
            }
        });
        openDialogHelpJButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                openDialogHelpJButtonActionPerformed(evt);
            }
        });

        advancedSettingsWarningLabel.setText("Click to open the Tide help page.");

        okButton.setText("OK");
        okButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                okButtonActionPerformed(evt);
            }
        });

        closeButton.setText("Close");
        closeButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                closeButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout backgroundPanelLayout = new javax.swing.GroupLayout(backgroundPanel);
        backgroundPanel.setLayout(backgroundPanelLayout);
        backgroundPanelLayout.setHorizontalGroup(
            backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, backgroundPanelLayout.createSequentialGroup()
                .addGap(20, 20, 20)
                .addComponent(openDialogHelpJButton)
                .addGap(18, 18, 18)
                .addComponent(advancedSettingsWarningLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(okButton, javax.swing.GroupLayout.PREFERRED_SIZE, 59, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(closeButton)
                .addGap(10, 10, 10))
            .addGroup(backgroundPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(tabbedPane, javax.swing.GroupLayout.PREFERRED_SIZE, 476, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        backgroundPanelLayout.setVerticalGroup(
            backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(backgroundPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(tabbedPane, javax.swing.GroupLayout.PREFERRED_SIZE, 426, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addGroup(backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.CENTER)
                    .addComponent(openDialogHelpJButton)
                    .addComponent(advancedSettingsWarningLabel)
                    .addComponent(okButton)
                    .addComponent(closeButton))
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
     * Close the dialog without saving the settings.
     *
     * @param evt
     */
    private void closeButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_closeButtonActionPerformed
        cancelled = true;
        dispose();
    }//GEN-LAST:event_closeButtonActionPerformed

    /**
     * Save the settings and then close the dialog.
     *
     * @param evt
     */
    private void okButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_okButtonActionPerformed
        if (validateInput(true)) {
            dispose();
        }
    }//GEN-LAST:event_okButtonActionPerformed

    /**
     * Change the cursor to a hand cursor.
     *
     * @param evt
     */
    private void openDialogHelpJButtonMouseEntered(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_openDialogHelpJButtonMouseEntered
        setCursor(new java.awt.Cursor(java.awt.Cursor.HAND_CURSOR));
    }//GEN-LAST:event_openDialogHelpJButtonMouseEntered

    /**
     * Change the cursor back to the default cursor.
     *
     * @param evt
     */
    private void openDialogHelpJButtonMouseExited(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_openDialogHelpJButtonMouseExited
        setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }//GEN-LAST:event_openDialogHelpJButtonMouseExited

    /**
     * Open the Tide help page.
     *
     * @param evt
     */
    private void openDialogHelpJButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_openDialogHelpJButtonActionPerformed
        setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        BareBonesBrowserLaunch.openURL("http://cruxtoolkit.sourceforge.net/");
        setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }//GEN-LAST:event_openDialogHelpJButtonActionPerformed

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void numberMatchesTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_numberMatchesTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_numberMatchesTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void maxPtmsPerPeptideTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_maxPtmsPerPeptideTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_maxPtmsPerPeptideTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void minSpectrumMzTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_minSpectrumMzTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_minSpectrumMzTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void maxSpectrumMzTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_maxSpectrumMzTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_maxSpectrumMzTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void mzBinWidthTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_mzBinWidthTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_mzBinWidthTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void minPeaksTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_minPeaksTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_minPeaksTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void removePrecursorPeakToleranceTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_removePrecursorPeakToleranceTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_removePrecursorPeakToleranceTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void minPepLengthTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_minPepLengthTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_minPepLengthTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void maxPepLengthTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_maxPepLengthTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_maxPepLengthTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void maxPrecursorMassTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_maxPrecursorMassTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_maxPrecursorMassTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void minPrecursorMassTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_minPrecursorMassTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_minPrecursorMassTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void decoySeedTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_decoySeedTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_decoySeedTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void maxVariablePtmsPerTypeTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_maxVariablePtmsPerTypeTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_maxVariablePtmsPerTypeTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void mzBinOffsetTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_mzBinOffsetTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_mzBinOffsetTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void mzBinWidthTxtActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_mzBinWidthTxtActionPerformed
        validateInput(false);
    }//GEN-LAST:event_mzBinWidthTxtActionPerformed

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void useFlankingCmbActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_useFlankingCmbActionPerformed
        validateInput(false);
    }//GEN-LAST:event_useFlankingCmbActionPerformed

    /**
     * Check if the output format is compatible with PeptideShaker.
     *
     * @param evt
     */
    private void outputFormatComboActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_outputFormatComboActionPerformed
        if (outputFormatCombo.getSelectedIndex() != 0 && this.isVisible()) {
            // invoke later to give time for components to update
            SwingUtilities.invokeLater(new Runnable() {
                public void run() {
                    JOptionPane.showMessageDialog(TideParametersDialog.this, JOptionEditorPane.getJOptionEditorPane(
                            "Note that the Tide " + (String) outputFormatCombo.getSelectedItem()
                            + " format is not compatible with <a href=\"https://compomics.github.io/projects/peptide-shaker.html\">PeptideShaker</a>."),
                            "Format Warning", JOptionPane.WARNING_MESSAGE);
                }
            });
        }
    }//GEN-LAST:event_outputFormatComboActionPerformed

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void minPtmsPerPeptideTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_minPtmsPerPeptideTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_minPtmsPerPeptideTxtKeyReleased

    /**
     * Inspects the parameter validity.
     *
     * @param showMessage if true an error messages are shown to the users
     * @return a boolean indicating if the parameters are valid
     */
    public boolean validateInput(boolean showMessage) {

        boolean valid = true;

        valid = GuiUtilities.validateIntegerInput(this, peptideLengthLabel, minPepLengthTxt, "minimum peptide length", "Peptide Length Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, peptideLengthLabel, maxPepLengthTxt, "minimum peptide length", "Peptide Length Error", true, showMessage, valid);
        valid = GuiUtilities.validateDoubleInput(this, precursorMassLabel, minPrecursorMassTxt, "minimum precursor mass", "Precursor Mass Error", true, showMessage, valid);
        valid = GuiUtilities.validateDoubleInput(this, precursorMassLabel, maxPrecursorMassTxt, "maximum precursor mass", "Precursor Mass Error", true, showMessage, valid);

        if (!minPtmsPerPeptideTxt.getText().trim().isEmpty()) {
            valid = GuiUtilities.validateIntegerInput(this, variablePtmsPerPeptideLabel, minPtmsPerPeptideTxt, "minimum number of variable PTMs", "Variable PTMs Error", true, showMessage, valid);
        }

        if (!maxPtmsPerPeptideTxt.getText().trim().isEmpty()) {
            valid = GuiUtilities.validateIntegerInput(this, variablePtmsPerPeptideLabel, maxPtmsPerPeptideTxt, "maximum number of variable PTMs", "Variable PTMs Error", true, showMessage, valid);
        }

        valid = GuiUtilities.validateIntegerInput(this, maxVariablePtmsPerTypeLabel, maxVariablePtmsPerTypeTxt, "maximum number of variable PTMs per type", "Variable PTMs Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, decoySeedLabel, decoySeedTxt, "decoy seed", "Decoy Seed Error", true, showMessage, valid);
        valid = GuiUtilities.validateDoubleInput(this, spectrumMzLabel, minSpectrumMzTxt, "minimum spectrum mz", "Spectrum Mz Error", true, showMessage, valid);
        if (!maxSpectrumMzTxt.getText().trim().isEmpty()) {
            valid = GuiUtilities.validateDoubleInput(this, spectrumMzLabel, maxSpectrumMzTxt, "maximum spectrum mz", "Spectrum Mz Error", true, showMessage, valid);
        }
        valid = GuiUtilities.validateIntegerInput(this, minPeaksLbl, minPeaksTxt, "minimum number of peaks", "Spectrum Peaks Error", true, showMessage, valid);
        valid = GuiUtilities.validateDoubleInput(this, removePrecursorPeakLabel, removePrecursorPeakToleranceTxt, "remove precursor peak tolerance", "Precursor Peak Error", true, showMessage, valid);
        valid = GuiUtilities.validateDoubleInput(this, mzBinWidthLabel, mzBinWidthTxt, "mz bin width", "Mz Bin Width Error", true, showMessage, valid);
        valid = GuiUtilities.validateDoubleInput(this, mzBinOffsetLabel, mzBinOffsetTxt, "mz bin offset", "Mz Bin Offset Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, numberMatchesLabel, numberMatchesTxt, "number of spectrum matches", "Number of Spectrum Matches Error", true, showMessage, valid);

        // peptide length: the low value should be lower than the high value
        try {
            double lowValue = Double.parseDouble(minPepLengthTxt.getText().trim());
            double highValue = Double.parseDouble(maxPepLengthTxt.getText().trim());

            if (lowValue > highValue) {
                if (showMessage && valid) {
                    JOptionPane.showMessageDialog(this, "The lower range value has to be smaller than the upper range value.",
                            "Peptide Length Error", JOptionPane.WARNING_MESSAGE);
                }
                valid = false;
                peptideLengthLabel.setForeground(Color.RED);
                peptideLengthLabel.setToolTipText("Please select a valid range (upper <= higher)");
            }
        } catch (NumberFormatException e) {
            // ignore, handled above
        }

        // precursor mass range: the low value should be lower than the high value
        try {
            double lowValue = Double.parseDouble(minPrecursorMassTxt.getText().trim());
            double highValue = Double.parseDouble(maxPrecursorMassTxt.getText().trim());

            if (lowValue > highValue) {
                if (showMessage && valid) {
                    JOptionPane.showMessageDialog(this, "The lower range value has to be smaller than the upper range value.",
                            "Precursor Mass Range Error", JOptionPane.WARNING_MESSAGE);
                }
                valid = false;
                precursorMassLabel.setForeground(Color.RED);
                precursorMassLabel.setToolTipText("Please select a valid range (upper <= higher)");
            }
        } catch (NumberFormatException e) {
            // ignore, handled above
        }

        // spectrum mz range: the low value should be lower than the high value
        try {
            double lowValue = Double.parseDouble(minSpectrumMzTxt.getText().trim());
            double highValue = Double.parseDouble(maxSpectrumMzTxt.getText().trim());

            if (lowValue > highValue) {
                if (showMessage && valid) {
                    JOptionPane.showMessageDialog(this, "The lower range value has to be smaller than the upper range value.",
                            "Spectrum Mz Range Error", JOptionPane.WARNING_MESSAGE);
                }
                valid = false;
                spectrumMzLabel.setForeground(Color.RED);
                spectrumMzLabel.setToolTipText("Please select a valid range (upper <= higher)");
            }
        } catch (NumberFormatException e) {
            // ignore, handled above
        }

        okButton.setEnabled(valid);
        return valid;
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel advancedSettingsWarningLabel;
    private javax.swing.JPanel backgroundPanel;
    private javax.swing.JButton closeButton;
    private javax.swing.JComboBox decoyFormatCombo;
    private javax.swing.JLabel decoySeedLabel;
    private javax.swing.JTextField decoySeedTxt;
    private javax.swing.JLabel decoyformatLabel;
    private javax.swing.JComboBox enzymeTypeCmb;
    private javax.swing.JLabel enzymeTypeLabel;
    private javax.swing.JComboBox exactPvalueCombo;
    private javax.swing.JLabel exactPvalueLabel;
    private javax.swing.JPanel indexPanel;
    private javax.swing.JComboBox keepTerminalAaCombo;
    private javax.swing.JLabel keepTerminalAaLabel;
    private javax.swing.JTextField maxPepLengthTxt;
    private javax.swing.JTextField maxPrecursorMassTxt;
    private javax.swing.JTextField maxPtmsPerPeptideTxt;
    private javax.swing.JTextField maxSpectrumMzTxt;
    private javax.swing.JLabel maxVariablePtmsPerTypeLabel;
    private javax.swing.JTextField maxVariablePtmsPerTypeTxt;
    private javax.swing.JLabel minPeaksLbl;
    private javax.swing.JTextField minPeaksTxt;
    private javax.swing.JTextField minPepLengthTxt;
    private javax.swing.JTextField minPrecursorMassTxt;
    private javax.swing.JTextField minPtmsPerPeptideTxt;
    private javax.swing.JTextField minSpectrumMzTxt;
    private javax.swing.JComboBox monoPrecursorCmb;
    private javax.swing.JLabel monoPrecursorLabel;
    private javax.swing.JLabel mzBinOffsetLabel;
    private javax.swing.JTextField mzBinOffsetTxt;
    private javax.swing.JLabel mzBinWidthLabel;
    private javax.swing.JTextField mzBinWidthTxt;
    private javax.swing.JLabel numberMatchesLabel;
    private javax.swing.JTextField numberMatchesTxt;
    private javax.swing.JButton okButton;
    private javax.swing.JButton openDialogHelpJButton;
    private javax.swing.JComboBox outputFormatCombo;
    private javax.swing.JLabel outputFormatLabel;
    private javax.swing.JLabel peptideLengthDividerLabel;
    private javax.swing.JLabel peptideLengthLabel;
    private javax.swing.JComboBox peptideListCmb;
    private javax.swing.JLabel peptideListLabel;
    private javax.swing.JLabel precursorMassDividerLabel;
    private javax.swing.JLabel precursorMassLabel;
    private javax.swing.JComboBox removeMethionineCmb;
    private javax.swing.JLabel removeMethionineLabel;
    private javax.swing.JComboBox removePrecursorPeakCombo;
    private javax.swing.JLabel removePrecursorPeakLabel;
    private javax.swing.JLabel removePrecursorPeakToleranceLbl;
    private javax.swing.JTextField removePrecursorPeakToleranceTxt;
    private javax.swing.JComboBox removeTempFoldersCmb;
    private javax.swing.JLabel removeTempFoldersLabel;
    private javax.swing.JComboBox scoreFunctionCombo;
    private javax.swing.JLabel scoreFunctionLabel;
    private javax.swing.JPanel searchPanel;
    private javax.swing.JLabel spectrumMzDividerLabel;
    private javax.swing.JLabel spectrumMzLabel;
    private javax.swing.JTabbedPane tabbedPane;
    private javax.swing.JComboBox useFlankingCmb;
    private javax.swing.JLabel useFlankingLabel;
    private javax.swing.JComboBox useNeutralLossCmb;
    private javax.swing.JLabel useNeutralLossLabel;
    private javax.swing.JLabel variablePtmsPerPeptideDividerLabel;
    private javax.swing.JLabel variablePtmsPerPeptideLabel;
    // End of variables declaration//GEN-END:variables
}
