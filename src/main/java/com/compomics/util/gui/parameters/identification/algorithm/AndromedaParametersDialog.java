package com.compomics.util.gui.parameters.identification.algorithm;

import com.compomics.util.examples.BareBonesBrowserLaunch;
import com.compomics.util.gui.parameters.identification.IdentificationAlgorithmParameter;
import com.compomics.util.parameters.identification.tool_specific.AndromedaParameters;
import com.compomics.util.parameters.identification.tool_specific.AndromedaParameters.AndromedaDecoyMode;
import com.compomics.util.experiment.mass_spectrometry.FragmentationMethod;
import com.compomics.util.gui.GuiUtilities;
import java.awt.Color;
import java.awt.Dialog;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JOptionPane;
import javax.swing.SwingConstants;
import com.compomics.util.gui.parameters.identification.AlgorithmParametersDialog;

/**
 * Dialog for the Andromeda specific settings.
 *
 * @author Harald Barsnes
 */
public class AndromedaParametersDialog extends javax.swing.JDialog implements AlgorithmParametersDialog {

    /**
     * Boolean indicating whether the used canceled the editing.
     */
    private boolean cancelled = false;
    /**
     * Boolean indicating whether the settings can be edited by the user.
     */
    private boolean editable;

    /**
     * Creates a new AndromedaParametersDialog with a frame as owner.
     *
     * @param parentFrame the parent frame
     * @param andromedaParameters the Andromeda parameters
     * @param editable boolean indicating whether the settings can be edited by
     * the user
     */
    public AndromedaParametersDialog(java.awt.Frame parentFrame, AndromedaParameters andromedaParameters, boolean editable) {
        super(parentFrame, true);
        this.editable = editable;
        initComponents();
        setUpGUI();
        populateGUI(andromedaParameters);
        validateInput(false);
        setLocationRelativeTo(parentFrame);
        setVisible(true);
    }

    /**
     * Creates a new AndromedaParametersDialog with a dialog as owner.
     *
     * @param owner the dialog owner
     * @param parentFrame the parent frame
     * @param andromedaParameters the Andromeda parameters
     * @param editable boolean indicating whether the settings can be edited by
     * the user
     */
    public AndromedaParametersDialog(Dialog owner, java.awt.Frame parentFrame, AndromedaParameters andromedaParameters, boolean editable) {
        super(owner, true);
        this.editable = editable;
        initComponents();
        setUpGUI();
        populateGUI(andromedaParameters);
        validateInput(false);
        setLocationRelativeTo(owner);
        setVisible(true);
    }

    /**
     * Sets up the GUI.
     */
    private void setUpGUI() {

        neutralLossesCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        fragMethodCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        fragmentAllCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        empiricalCorrectionCombo.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        higherChargeCombo.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        ammoniaLossCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        waterLossCombo.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        equalILCombo.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        decoyModeCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));

        minPepLengthNoEnzymeTxt.setEditable(editable);
        minPepLengthNoEnzymeTxt.setEnabled(editable);
        maxPepLengthNoEnzymeTxt.setEditable(editable);
        maxPepLengthNoEnzymeTxt.setEnabled(editable);
        maxPeptideMassTxt.setEditable(editable);
        maxPeptideMassTxt.setEnabled(editable);
        numberMatchesTxt.setEditable(editable);
        numberMatchesTxt.setEnabled(editable);
        maxPtmsTxt.setEditable(editable);
        maxPtmsTxt.setEnabled(editable);
        fragMethodCmb.setEnabled(editable);
        waterLossCombo.setEnabled(editable);
        ammoniaLossCmb.setEnabled(editable);
        neutralLossesCmb.setEnabled(editable);
        equalILCombo.setEnabled(editable);
        fragmentAllCmb.setEnabled(editable);
        empiricalCorrectionCombo.setEnabled(editable);
        higherChargeCombo.setEnabled(editable);
        maxCombinationsTxt.setEditable(editable);
        maxCombinationsTxt.setEnabled(editable);
        topPeaksTxt.setEditable(editable);
        topPeaksTxt.setEnabled(editable);
        topPeaksWindowTxt.setEditable(editable);
        topPeaksWindowTxt.setEnabled(editable);
        decoyModeCmb.setEnabled(editable);
    }

    /**
     * Populates the GUI using the given settings.
     *
     * @param andromedaParameters the parameters to display
     */
    private void populateGUI(AndromedaParameters andromedaParameters) {

        minPepLengthNoEnzymeTxt.setText(andromedaParameters.getMinPeptideLengthNoEnzyme() + "");
        maxPepLengthNoEnzymeTxt.setText(andromedaParameters.getMaxPeptideLengthNoEnzyme() + "");
        maxPeptideMassTxt.setText(andromedaParameters.getMaxPeptideMass() + "");
        numberMatchesTxt.setText(andromedaParameters.getNumberOfCandidates() + "");
        maxPtmsTxt.setText(andromedaParameters.getMaxNumberOfModifications() + "");

        if (andromedaParameters.getFragmentationMethod() == FragmentationMethod.CID) {
            fragMethodCmb.setSelectedIndex(0);
        } else if (andromedaParameters.getFragmentationMethod() == FragmentationMethod.HCD) {
            fragMethodCmb.setSelectedIndex(1);
        } else if (andromedaParameters.getFragmentationMethod() == FragmentationMethod.ETD) {
            fragMethodCmb.setSelectedIndex(2);
        }

        if (andromedaParameters.isIncludeWater()) {
            waterLossCombo.setSelectedIndex(0);
        } else {
            waterLossCombo.setSelectedIndex(1);
        }
        if (andromedaParameters.isIncludeAmmonia()) {
            ammoniaLossCmb.setSelectedIndex(0);
        } else {
            ammoniaLossCmb.setSelectedIndex(1);
        }
        if (andromedaParameters.isDependentLosses()) {
            neutralLossesCmb.setSelectedIndex(0);
        } else {
            neutralLossesCmb.setSelectedIndex(1);
        }
        if (andromedaParameters.isEqualIL()) {
            equalILCombo.setSelectedIndex(0);
        } else {
            equalILCombo.setSelectedIndex(1);
        }
        if (andromedaParameters.isFragmentAll()) {
            fragmentAllCmb.setSelectedIndex(0);
        } else {
            fragmentAllCmb.setSelectedIndex(1);
        }
        if (andromedaParameters.isEmpiricalCorrection()) {
            empiricalCorrectionCombo.setSelectedIndex(0);
        } else {
            empiricalCorrectionCombo.setSelectedIndex(1);
        }
        if (andromedaParameters.isHigherCharge()) {
            higherChargeCombo.setSelectedIndex(0);
        } else {
            higherChargeCombo.setSelectedIndex(1);
        }

        decoyModeCmb.setSelectedItem(andromedaParameters.getDecoyMode());

        maxCombinationsTxt.setText(andromedaParameters.getMaxCombinations() + "");
        topPeaksTxt.setText(andromedaParameters.getTopPeaks() + "");
        topPeaksWindowTxt.setText(andromedaParameters.getTopPeaksWindow() + "");
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
     * Returns the user selection as Andromeda parameters object.
     *
     * @return the user selection
     */
    public AndromedaParameters getInput() {

        AndromedaParameters result = new AndromedaParameters();

        String input = minPepLengthNoEnzymeTxt.getText().trim();
        if (!input.equals("")) {
            result.setMinPeptideLengthNoEnzyme(new Integer(input));
        }
        input = maxPepLengthNoEnzymeTxt.getText().trim();
        if (!input.equals("")) {
            result.setMaxPeptideLengthNoEnzyme(new Integer(input));
        }
        input = maxPeptideMassTxt.getText().trim();
        if (!input.equals("")) {
            result.setMaxPeptideMass(new Double(input));
        }
        input = numberMatchesTxt.getText().trim();
        if (!input.equals("")) {
            result.setNumberOfCandidates(new Integer(input));
        }
        input = maxPtmsTxt.getText().trim();
        if (!input.equals("")) {
            result.setMaxNumberOfModifications(new Integer(input));
        }

        switch (fragMethodCmb.getSelectedIndex()) {
            case 0:
                result.setFragmentationMethod(FragmentationMethod.CID);
                break;
            case 1:
                result.setFragmentationMethod(FragmentationMethod.HCD);
                break;
            default:
                result.setFragmentationMethod(FragmentationMethod.ETD);
                break;
        }

        result.setIncludeWater(waterLossCombo.getSelectedIndex() == 0);
        result.setIncludeAmmonia(ammoniaLossCmb.getSelectedIndex() == 0);
        result.setDependentLosses(neutralLossesCmb.getSelectedIndex() == 0);
        result.setEqualIL(equalILCombo.getSelectedIndex() == 0);
        result.setFragmentAll(fragmentAllCmb.getSelectedIndex() == 0);
        result.setEmpiricalCorrection(empiricalCorrectionCombo.getSelectedIndex() == 0);
        result.setHigherCharge(higherChargeCombo.getSelectedIndex() == 0);

        input = maxCombinationsTxt.getText().trim();
        if (!input.equals("")) {
            result.setMaxCombinations(new Integer(input));
        }
        input = topPeaksTxt.getText().trim();
        if (!input.equals("")) {
            result.setTopPeaks(new Integer(input));
        }
        input = topPeaksWindowTxt.getText().trim();
        if (!input.equals("")) {
            result.setTopPeaksWindow(new Integer(input));
        }

        result.setDecoyMode((AndromedaDecoyMode) decoyModeCmb.getSelectedItem());

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
        settingsPanel = new javax.swing.JPanel();
        settingsScrollPane = new javax.swing.JScrollPane();
        settingsInnerPanel = new javax.swing.JPanel();
        peptideLengthNoEnzymeLabel = new javax.swing.JLabel();
        minPepLengthNoEnzymeTxt = new javax.swing.JTextField();
        peptideLengthNoEnzymeDividerLabel = new javax.swing.JLabel();
        maxPepLengthNoEnzymeTxt = new javax.swing.JTextField();
        maxPeptideMassLabel = new javax.swing.JLabel();
        maxPeptideMassTxt = new javax.swing.JTextField();
        numberMatchesLabel = new javax.swing.JLabel();
        numberMatchesTxt = new javax.swing.JTextField();
        fragMethodLabel = new javax.swing.JLabel();
        fragMethodCmb = new javax.swing.JComboBox();
        ammoniaLossLabel = new javax.swing.JLabel();
        ammoniaLossCmb = new javax.swing.JComboBox();
        maxPtmsLabel = new javax.swing.JLabel();
        maxPtmsTxt = new javax.swing.JTextField();
        empiricalCorrectionLabel = new javax.swing.JLabel();
        empiricalCorrectionCombo = new javax.swing.JComboBox();
        higherChargeLabel = new javax.swing.JLabel();
        higherChargeCombo = new javax.swing.JComboBox();
        topPeaksLabel = new javax.swing.JLabel();
        topPeaksTxt = new javax.swing.JTextField();
        neutralLossesLabel = new javax.swing.JLabel();
        neutralLossesCmb = new javax.swing.JComboBox();
        fragmentAllLabel = new javax.swing.JLabel();
        fragmentAllCmb = new javax.swing.JComboBox();
        maxCombinationsLabel = new javax.swing.JLabel();
        maxCombinationsTxt = new javax.swing.JTextField();
        topPeaksWindowLabel = new javax.swing.JLabel();
        topPeaksWindowTxt = new javax.swing.JTextField();
        waterLossLabel = new javax.swing.JLabel();
        waterLossCombo = new javax.swing.JComboBox();
        equalILLabel = new javax.swing.JLabel();
        equalILCombo = new javax.swing.JComboBox();
        decoyModeLabel = new javax.swing.JLabel();
        decoyModeCmb = new javax.swing.JComboBox();
        openDialogHelpJButton = new javax.swing.JButton();
        advancedSettingsWarningLabel = new javax.swing.JLabel();
        okButton = new javax.swing.JButton();
        closeButton = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Andromeda Advanced Settings");
        setMinimumSize(new java.awt.Dimension(400, 400));

        backgroundPanel.setBackground(new java.awt.Color(230, 230, 230));

        settingsPanel.setBackground(new java.awt.Color(230, 230, 230));
        settingsPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Search Settings"));
        settingsPanel.setPreferredSize(new java.awt.Dimension(518, 143));

        settingsScrollPane.setBorder(null);

        settingsInnerPanel.setBackground(new java.awt.Color(230, 230, 230));

        peptideLengthNoEnzymeLabel.setText("Peptide Length No Enzyme");

        minPepLengthNoEnzymeTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        minPepLengthNoEnzymeTxt.setText("8");
        minPepLengthNoEnzymeTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                minPepLengthNoEnzymeTxtKeyReleased(evt);
            }
        });

        peptideLengthNoEnzymeDividerLabel.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        peptideLengthNoEnzymeDividerLabel.setText("-");

        maxPepLengthNoEnzymeTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        maxPepLengthNoEnzymeTxt.setText("30");
        maxPepLengthNoEnzymeTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                maxPepLengthNoEnzymeTxtKeyReleased(evt);
            }
        });

        maxPeptideMassLabel.setText("Max Peptide Mass");

        maxPeptideMassTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        maxPeptideMassTxt.setText("4600");
        maxPeptideMassTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                maxPeptideMassTxtKeyReleased(evt);
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

        fragMethodLabel.setText("Fragmentation Method");

        fragMethodCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "CID", "HCD", "ETD" }));

        ammoniaLossLabel.setText("Ammonia Loss");

        ammoniaLossCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));

        maxPtmsLabel.setText("Max Variable PTMs");

        maxPtmsTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        maxPtmsTxt.setText("5");
        maxPtmsTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                maxPtmsTxtKeyReleased(evt);
            }
        });

        empiricalCorrectionLabel.setText("Empirical Correction");

        empiricalCorrectionCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));

        higherChargeLabel.setText("Higher Charge");

        higherChargeCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));

        topPeaksLabel.setText("Top Peaks");

        topPeaksTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        topPeaksTxt.setText("8");
        topPeaksTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                topPeaksTxtKeyReleased(evt);
            }
        });

        neutralLossesLabel.setText("Sequence Dependent Neutral Loss");

        neutralLossesCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));

        fragmentAllLabel.setText("Fragment All");

        fragmentAllCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));

        maxCombinationsLabel.setText("Max Combinations");

        maxCombinationsTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        maxCombinationsTxt.setText("250");
        maxCombinationsTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                maxCombinationsTxtKeyReleased(evt);
            }
        });

        topPeaksWindowLabel.setText("Top Peaks Windows");

        topPeaksWindowTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        topPeaksWindowTxt.setText("100");
        topPeaksWindowTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                topPeaksWindowTxtKeyReleased(evt);
            }
        });

        waterLossLabel.setText("Water Loss");

        waterLossCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));

        equalILLabel.setText("Equal I and L");

        equalILCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));
        equalILCombo.setSelectedIndex(1);

        decoyModeLabel.setText("Decoy Mode");

        decoyModeCmb.setModel(new DefaultComboBoxModel(AndromedaDecoyMode.values()));

        javax.swing.GroupLayout settingsInnerPanelLayout = new javax.swing.GroupLayout(settingsInnerPanel);
        settingsInnerPanel.setLayout(settingsInnerPanelLayout);
        settingsInnerPanelLayout.setHorizontalGroup(
            settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(settingsInnerPanelLayout.createSequentialGroup()
                .addComponent(peptideLengthNoEnzymeLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(21, 21, 21)
                .addComponent(minPepLengthNoEnzymeTxt, javax.swing.GroupLayout.DEFAULT_SIZE, 78, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(peptideLengthNoEnzymeDividerLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 27, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(7, 7, 7)
                .addComponent(maxPepLengthNoEnzymeTxt, javax.swing.GroupLayout.DEFAULT_SIZE, 84, Short.MAX_VALUE))
            .addGroup(settingsInnerPanelLayout.createSequentialGroup()
                .addComponent(maxPeptideMassLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(21, 21, 21)
                .addComponent(maxPeptideMassTxt))
            .addGroup(settingsInnerPanelLayout.createSequentialGroup()
                .addComponent(numberMatchesLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(21, 21, 21)
                .addComponent(numberMatchesTxt))
            .addGroup(settingsInnerPanelLayout.createSequentialGroup()
                .addComponent(neutralLossesLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(21, 21, 21)
                .addComponent(neutralLossesCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            .addGroup(settingsInnerPanelLayout.createSequentialGroup()
                .addComponent(ammoniaLossLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(21, 21, 21)
                .addComponent(ammoniaLossCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            .addGroup(settingsInnerPanelLayout.createSequentialGroup()
                .addComponent(waterLossLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(21, 21, 21)
                .addComponent(waterLossCombo, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, settingsInnerPanelLayout.createSequentialGroup()
                .addComponent(equalILLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(21, 21, 21)
                .addComponent(equalILCombo, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            .addGroup(settingsInnerPanelLayout.createSequentialGroup()
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(maxPtmsLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(fragMethodLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(21, 21, 21)
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(fragMethodCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(maxPtmsTxt)))
            .addGroup(settingsInnerPanelLayout.createSequentialGroup()
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(fragmentAllLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(empiricalCorrectionLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(higherChargeLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(21, 21, 21)
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(empiricalCorrectionCombo, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(higherChargeCombo, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(fragmentAllCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
            .addGroup(settingsInnerPanelLayout.createSequentialGroup()
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                    .addComponent(topPeaksWindowLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(maxCombinationsLabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(topPeaksLabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(21, 21, 21)
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(topPeaksTxt)
                    .addComponent(topPeaksWindowTxt)
                    .addComponent(maxCombinationsTxt)))
            .addGroup(settingsInnerPanelLayout.createSequentialGroup()
                .addComponent(decoyModeLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(21, 21, 21)
                .addComponent(decoyModeCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        settingsInnerPanelLayout.setVerticalGroup(
            settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(settingsInnerPanelLayout.createSequentialGroup()
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(minPepLengthNoEnzymeTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(maxPepLengthNoEnzymeTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(peptideLengthNoEnzymeDividerLabel)
                    .addComponent(peptideLengthNoEnzymeLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(maxPeptideMassTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(maxPeptideMassLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(numberMatchesLabel)
                    .addComponent(numberMatchesTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(maxPtmsLabel)
                    .addComponent(maxPtmsTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fragMethodLabel)
                    .addComponent(fragMethodCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(waterLossLabel)
                    .addComponent(waterLossCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(ammoniaLossLabel)
                    .addComponent(ammoniaLossCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(neutralLossesLabel)
                    .addComponent(neutralLossesCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(equalILCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(equalILLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fragmentAllLabel)
                    .addComponent(fragmentAllCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(empiricalCorrectionLabel)
                    .addComponent(empiricalCorrectionCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(higherChargeLabel)
                    .addComponent(higherChargeCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(maxCombinationsLabel)
                    .addComponent(maxCombinationsTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(topPeaksLabel)
                    .addComponent(topPeaksTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(topPeaksWindowTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(topPeaksWindowLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(settingsInnerPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(decoyModeLabel)
                    .addComponent(decoyModeCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap())
        );

        settingsScrollPane.setViewportView(settingsInnerPanel);

        javax.swing.GroupLayout settingsPanelLayout = new javax.swing.GroupLayout(settingsPanel);
        settingsPanel.setLayout(settingsPanelLayout);
        settingsPanelLayout.setHorizontalGroup(
            settingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(settingsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(settingsScrollPane)
                .addContainerGap())
        );
        settingsPanelLayout.setVerticalGroup(
            settingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(settingsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(settingsScrollPane, javax.swing.GroupLayout.DEFAULT_SIZE, 540, Short.MAX_VALUE)
                .addContainerGap())
        );

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

        advancedSettingsWarningLabel.setText("Click to open the Andromeda help page.");

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
            .addGroup(backgroundPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(settingsPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 481, Short.MAX_VALUE)
                    .addGroup(backgroundPanelLayout.createSequentialGroup()
                        .addGap(10, 10, 10)
                        .addComponent(openDialogHelpJButton)
                        .addGap(18, 18, 18)
                        .addComponent(advancedSettingsWarningLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 92, Short.MAX_VALUE)
                        .addComponent(okButton, javax.swing.GroupLayout.PREFERRED_SIZE, 59, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(closeButton)))
                .addContainerGap())
        );
        backgroundPanelLayout.setVerticalGroup(
            backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(backgroundPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(settingsPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 585, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
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
     * Open the Andromeda help page.
     *
     * @param evt
     */
    private void openDialogHelpJButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_openDialogHelpJButtonActionPerformed
        setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        BareBonesBrowserLaunch.openURL("http://coxdocs.org/doku.php?id=maxquant:andromeda:start");
        setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }//GEN-LAST:event_openDialogHelpJButtonActionPerformed

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void topPeaksWindowTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_topPeaksWindowTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_topPeaksWindowTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void maxPtmsTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_maxPtmsTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_maxPtmsTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void minPepLengthNoEnzymeTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_minPepLengthNoEnzymeTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_minPepLengthNoEnzymeTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void maxPepLengthNoEnzymeTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_maxPepLengthNoEnzymeTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_maxPepLengthNoEnzymeTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void maxPeptideMassTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_maxPeptideMassTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_maxPeptideMassTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void topPeaksTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_topPeaksTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_topPeaksTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void maxCombinationsTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_maxCombinationsTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_maxCombinationsTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void numberMatchesTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_numberMatchesTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_numberMatchesTxtKeyReleased

    /**
     * Inspects the parameter validity.
     *
     * @param showMessage if true an error messages are shown to the users
     * @return a boolean indicating if the parameters are valid
     */
    public boolean validateInput(boolean showMessage) {

        boolean valid = true;

        valid = GuiUtilities.validateIntegerInput(this, peptideLengthNoEnzymeLabel, minPepLengthNoEnzymeTxt, "minimum peptide length", "Peptide Length Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, peptideLengthNoEnzymeLabel, maxPepLengthNoEnzymeTxt, "maximum peptide length", "Peptide Length Error", true, showMessage, valid);
        valid = GuiUtilities.validateDoubleInput(this, maxPeptideMassLabel, maxPeptideMassTxt, "maximum peptide mass", "Peptide Mass Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, numberMatchesLabel, numberMatchesTxt, "number of spectrum matches", "Number of Spectrum Matches Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, maxPtmsLabel, maxPtmsTxt, "maximum number of variable PTMs", "Variable PTMs Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, maxCombinationsLabel, maxCombinationsTxt, "maximum combinations", "Max Combinations Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, topPeaksLabel, topPeaksTxt, "top peaks", "Top Peaks Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, topPeaksWindowLabel, topPeaksWindowTxt, "top peaks window", "Top Peaks Window Error", true, showMessage, valid);

        // peptide length: the low value should be lower than the high value
        try {
            double lowValue = Double.parseDouble(minPepLengthNoEnzymeTxt.getText().trim());
            double highValue = Double.parseDouble(maxPepLengthNoEnzymeTxt.getText().trim());

            if (lowValue > highValue) {
                if (showMessage && valid) {
                    JOptionPane.showMessageDialog(this, "The lower range value has to be smaller than the upper range value.",
                            "Peptide Length Error", JOptionPane.WARNING_MESSAGE);
                }
                valid = false;
                peptideLengthNoEnzymeLabel.setForeground(Color.RED);
                peptideLengthNoEnzymeLabel.setToolTipText("Please select a valid range (upper <= higher)");
            }
        } catch (NumberFormatException e) {
            // ignore, handled above
        }

        okButton.setEnabled(valid);
        return valid;
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel advancedSettingsWarningLabel;
    private javax.swing.JComboBox ammoniaLossCmb;
    private javax.swing.JLabel ammoniaLossLabel;
    private javax.swing.JPanel backgroundPanel;
    private javax.swing.JButton closeButton;
    private javax.swing.JComboBox decoyModeCmb;
    private javax.swing.JLabel decoyModeLabel;
    private javax.swing.JComboBox empiricalCorrectionCombo;
    private javax.swing.JLabel empiricalCorrectionLabel;
    private javax.swing.JComboBox equalILCombo;
    private javax.swing.JLabel equalILLabel;
    private javax.swing.JComboBox fragMethodCmb;
    private javax.swing.JLabel fragMethodLabel;
    private javax.swing.JComboBox fragmentAllCmb;
    private javax.swing.JLabel fragmentAllLabel;
    private javax.swing.JComboBox higherChargeCombo;
    private javax.swing.JLabel higherChargeLabel;
    private javax.swing.JLabel maxCombinationsLabel;
    private javax.swing.JTextField maxCombinationsTxt;
    private javax.swing.JTextField maxPepLengthNoEnzymeTxt;
    private javax.swing.JLabel maxPeptideMassLabel;
    private javax.swing.JTextField maxPeptideMassTxt;
    private javax.swing.JLabel maxPtmsLabel;
    private javax.swing.JTextField maxPtmsTxt;
    private javax.swing.JTextField minPepLengthNoEnzymeTxt;
    private javax.swing.JComboBox neutralLossesCmb;
    private javax.swing.JLabel neutralLossesLabel;
    private javax.swing.JLabel numberMatchesLabel;
    private javax.swing.JTextField numberMatchesTxt;
    private javax.swing.JButton okButton;
    private javax.swing.JButton openDialogHelpJButton;
    private javax.swing.JLabel peptideLengthNoEnzymeDividerLabel;
    private javax.swing.JLabel peptideLengthNoEnzymeLabel;
    private javax.swing.JPanel settingsInnerPanel;
    private javax.swing.JPanel settingsPanel;
    private javax.swing.JScrollPane settingsScrollPane;
    private javax.swing.JLabel topPeaksLabel;
    private javax.swing.JTextField topPeaksTxt;
    private javax.swing.JLabel topPeaksWindowLabel;
    private javax.swing.JTextField topPeaksWindowTxt;
    private javax.swing.JComboBox waterLossCombo;
    private javax.swing.JLabel waterLossLabel;
    // End of variables declaration//GEN-END:variables

}
