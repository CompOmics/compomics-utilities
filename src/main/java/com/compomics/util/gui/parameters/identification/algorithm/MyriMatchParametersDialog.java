package com.compomics.util.gui.parameters.identification.algorithm;

import com.compomics.util.examples.BareBonesBrowserLaunch;
import com.compomics.util.gui.parameters.identification.IdentificationAlgorithmParameter;
import com.compomics.util.parameters.identification.tool_specific.MyriMatchParameters;
import com.compomics.util.gui.GuiUtilities;
import com.compomics.util.gui.JOptionEditorPane;
import java.awt.Color;
import java.awt.Dialog;
import javax.swing.JOptionPane;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import com.compomics.util.gui.parameters.identification.AlgorithmParametersDialog;

/**
 * Dialog for the MyriMatch specific settings.
 *
 * @author Harald Barsnes
 */
public class MyriMatchParametersDialog extends javax.swing.JDialog implements AlgorithmParametersDialog {

    /**
     * Empty default constructor
     */
    public MyriMatchParametersDialog() {
    }

    /**
     * Boolean indicating whether the used canceled the editing.
     */
    private boolean cancelled = false;
    /**
     * Boolean indicating whether the settings can be edited by the user.
     */
    private boolean editable;

    /**
     * Creates new form MyriMatchSettingsDialog with a frame as owner.
     *
     * @param parent the parent frame
     * @param myriMatchParameters the MyriMatch parameters
     * @param editable boolean indicating whether the settings can be edited by the user
     */
    public MyriMatchParametersDialog(java.awt.Frame parent, MyriMatchParameters myriMatchParameters, boolean editable) {
        super(parent, true);
        this.editable = editable;
        initComponents();
        setUpGUI();
        populateGUI(myriMatchParameters);
        validateInput(false);
        setLocationRelativeTo(parent);
        setVisible(true);
    }

    /**
     * Creates new form MyriMatchSettingsDialog with a dialog as owner.
     *
     * @param owner the dialog owner
     * @param parent the parent frame
     * @param myriMatchParameters the MyriMatch parameters
     * @param editable boolean indicating whether the settings can be edited by the user
     */
    public MyriMatchParametersDialog(Dialog owner, java.awt.Frame parent, MyriMatchParameters myriMatchParameters, boolean editable) {
        super(owner, true);
        this.editable = editable;
        initComponents();
        setUpGUI();
        populateGUI(myriMatchParameters);
        validateInput(false);
        setLocationRelativeTo(owner);
        setVisible(true);
    }

    /**
     * Sets up the GUI.
     */
    private void setUpGUI() {
        useSmartPlusThreeModelCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        computeXCorrCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        terminiCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        fragmentationMethodCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        outputFormatCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
    }

    /**
     * Populates the GUI using the given settings.
     * 
     * @param myriMatchParameters the parameters to display
     */
    private void populateGUI(MyriMatchParameters myriMatchParameters) {

        if (myriMatchParameters.getMinPeptideLength() != null) {
            minPepLengthTxt.setText(myriMatchParameters.getMinPeptideLength() + "");
        }
        if (myriMatchParameters.getMaxPeptideLength() != null) {
            maxPepLengthTxt.setText(myriMatchParameters.getMaxPeptideLength() + "");
        }

        if (myriMatchParameters.getMinPrecursorMass() != null) {
            minPrecursorMassTxt.setText(myriMatchParameters.getMinPrecursorMass() + "");
        }
        if (myriMatchParameters.getMaxPrecursorMass() != null) {
            maxPrecursorMassTxt.setText(myriMatchParameters.getMaxPrecursorMass() + "");
        }

        if (myriMatchParameters.getNumberOfSpectrumMatches() != null) {
            numberMatchesTxt.setText(myriMatchParameters.getNumberOfSpectrumMatches() + "");
        }

        if (myriMatchParameters.getMaxDynamicMods() != null) {
            maxPtmsTxt.setText(myriMatchParameters.getMaxDynamicMods() + "");
        }

        fragmentationMethodCmb.setSelectedItem(myriMatchParameters.getFragmentationRule());

        if (myriMatchParameters.getMinTerminiCleavages() != null) {
            terminiCmb.setSelectedIndex(myriMatchParameters.getMinTerminiCleavages());
        }

        if (myriMatchParameters.getUseSmartPlusThreeModel()) {
            useSmartPlusThreeModelCmb.setSelectedIndex(0);
        } else {
            useSmartPlusThreeModelCmb.setSelectedIndex(1);
        }

        if (myriMatchParameters.getComputeXCorr()) {
            computeXCorrCmb.setSelectedIndex(0);
        } else {
            computeXCorrCmb.setSelectedIndex(1);
        }

        if (myriMatchParameters.getTicCutoffPercentage() != null) {
            ticCutoffPercentageTxt.setText(myriMatchParameters.getTicCutoffPercentage() + "");
        }

        if (myriMatchParameters.getNumIntensityClasses() != null) {
            numIntensityClassesTxt.setText(myriMatchParameters.getNumIntensityClasses() + "");
        }

        if (myriMatchParameters.getClassSizeMultiplier() != null) {
            classSizeMultiplierTxt.setText(myriMatchParameters.getClassSizeMultiplier() + "");
        }

        if (myriMatchParameters.getNumberOfBatches() != null) {
            numbBatchesTxt.setText(myriMatchParameters.getNumberOfBatches() + "");
        }

        if (myriMatchParameters.getMaxPeakCount() != null) {
            maxPeakCountTxt.setText(myriMatchParameters.getMaxPeakCount() + "");
        }

        if (myriMatchParameters.getOutputFormat().equalsIgnoreCase("mzIdentML")) {
            outputFormatCmb.setSelectedIndex(0);
        } else {
            outputFormatCmb.setSelectedIndex(1);
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
     * Returns the user selection as MyriMatch parameters object.
     *
     * @return the user selection
     */
    public MyriMatchParameters getInput() {

        MyriMatchParameters result = new MyriMatchParameters();

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

        input = numberMatchesTxt.getText().trim();
        if (!input.equals("")) {
            result.setNumberOfSpectrumMatches(Integer.valueOf(input));
        }

        input = maxPtmsTxt.getText().trim();
        if (!input.equals("")) {
            result.setMaxDynamicMods(Integer.valueOf(input));
        }

        result.setFragmentationRule((String) fragmentationMethodCmb.getSelectedItem());
        result.setMinTerminiCleavages(terminiCmb.getSelectedIndex());
        result.setUseSmartPlusThreeModel(useSmartPlusThreeModelCmb.getSelectedIndex() == 0);
        result.setComputeXCorr(computeXCorrCmb.getSelectedIndex() == 0);

        input = ticCutoffPercentageTxt.getText().trim();
        if (!input.equals("")) {
            result.setTicCutoffPercentage(Double.valueOf(input));
        }

        input = numIntensityClassesTxt.getText().trim();
        if (!input.equals("")) {
            result.setNumIntensityClasses(Integer.valueOf(input));
        }

        input = classSizeMultiplierTxt.getText().trim();
        if (!input.equals("")) {
            result.setClassSizeMultiplier(Integer.valueOf(input));
        }

        input = numbBatchesTxt.getText().trim();
        if (!input.equals("")) {
            result.setNumberOfBatches(Integer.valueOf(input));
        }

        input = maxPeakCountTxt.getText().trim();
        if (!input.equals("")) {
            result.setMaxPeakCount(Integer.valueOf(input));
        }

        result.setOutputFormat((String) outputFormatCmb.getSelectedItem());

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
        advancedSearchSettingsPanel = new javax.swing.JPanel();
        scrollPane = new javax.swing.JScrollPane();
        advancedSettingsPanel = new javax.swing.JPanel();
        useSmartPlusThreeModelCmb = new javax.swing.JComboBox();
        useSmartPlusThreeModelLabel = new javax.swing.JLabel();
        minPepLengthTxt = new javax.swing.JTextField();
        peptideLengthDividerLabel = new javax.swing.JLabel();
        maxPepLengthTxt = new javax.swing.JTextField();
        peptideLengthLabel = new javax.swing.JLabel();
        numberMatchesLabel = new javax.swing.JLabel();
        numberMatchesTxt = new javax.swing.JTextField();
        computeXCorrlLabel = new javax.swing.JLabel();
        computeXCorrCmb = new javax.swing.JComboBox();
        numberTerminiLabel = new javax.swing.JLabel();
        maxPtmsLabel = new javax.swing.JLabel();
        maxPtmsTxt = new javax.swing.JTextField();
        terminiCmb = new javax.swing.JComboBox();
        precursorMassLabel = new javax.swing.JLabel();
        minPrecursorMassTxt = new javax.swing.JTextField();
        precursorMassDividerLabel = new javax.swing.JLabel();
        maxPrecursorMassTxt = new javax.swing.JTextField();
        ticCutoffPercentageLabel = new javax.swing.JLabel();
        ticCutoffPercentageTxt = new javax.swing.JTextField();
        numIntensityClassesLabel = new javax.swing.JLabel();
        numIntensityClassesTxt = new javax.swing.JTextField();
        classSizeMultiplierLabel = new javax.swing.JLabel();
        classSizeMultiplierTxt = new javax.swing.JTextField();
        numbBatchesLabel = new javax.swing.JLabel();
        numbBatchesTxt = new javax.swing.JTextField();
        fragmentationMethodLabel = new javax.swing.JLabel();
        fragmentationMethodCmb = new javax.swing.JComboBox();
        maxPeakCountLabel = new javax.swing.JLabel();
        maxPeakCountTxt = new javax.swing.JTextField();
        outputFormatLabel = new javax.swing.JLabel();
        outputFormatCmb = new javax.swing.JComboBox();
        okButton = new javax.swing.JButton();
        closeButton = new javax.swing.JButton();
        openDialogHelpJButton = new javax.swing.JButton();
        advancedSettingsWarningLabel = new javax.swing.JLabel();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("MyriMatch Advanced Settings");
        setMinimumSize(new java.awt.Dimension(400, 400));

        backgroundPanel.setBackground(new java.awt.Color(230, 230, 230));

        advancedSearchSettingsPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Search Settings"));
        advancedSearchSettingsPanel.setOpaque(false);

        scrollPane.setBorder(null);

        advancedSettingsPanel.setBackground(new java.awt.Color(230, 230, 230));

        useSmartPlusThreeModelCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));

        useSmartPlusThreeModelLabel.setText("Use Smart Plus Thee Model");

        minPepLengthTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        minPepLengthTxt.setText("8");
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

        peptideLengthLabel.setText("Peptide Length (min - max)");

        numberMatchesLabel.setText("Number of Spectrum Matches");

        numberMatchesTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        numberMatchesTxt.setText("1");
        numberMatchesTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                numberMatchesTxtKeyReleased(evt);
            }
        });

        computeXCorrlLabel.setText("Compute XCorr");

        computeXCorrCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));
        computeXCorrCmb.setSelectedIndex(1);

        numberTerminiLabel.setText("Enzymatic Terminals");

        maxPtmsLabel.setText("Max Variable PTMs per Peptide");

        maxPtmsTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        maxPtmsTxt.setText("2");
        maxPtmsTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                maxPtmsTxtKeyReleased(evt);
            }
        });

        terminiCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "None Required", "At Least One", "Both" }));
        terminiCmb.setSelectedIndex(2);

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

        ticCutoffPercentageLabel.setText("TIC Cutoff Percentage (0.0-1.0)");

        ticCutoffPercentageTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        ticCutoffPercentageTxt.setText("0.98");
        ticCutoffPercentageTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                ticCutoffPercentageTxtKeyReleased(evt);
            }
        });

        numIntensityClassesLabel.setText("Number of Intensity Classes");

        numIntensityClassesTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        numIntensityClassesTxt.setText("3");
        numIntensityClassesTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                numIntensityClassesTxtKeyReleased(evt);
            }
        });

        classSizeMultiplierLabel.setText("Class Size Multiplier");

        classSizeMultiplierTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        classSizeMultiplierTxt.setText("2");
        classSizeMultiplierTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                classSizeMultiplierTxtKeyReleased(evt);
            }
        });

        numbBatchesLabel.setText("Number of Batches");

        numbBatchesTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        numbBatchesTxt.setText("50");
        numbBatchesTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                numbBatchesTxtKeyReleased(evt);
            }
        });

        fragmentationMethodLabel.setText("Fragmentation Method");

        fragmentationMethodCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "CID", "HCD", "ETD" }));

        maxPeakCountLabel.setText("Max Peak Count");

        maxPeakCountTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        maxPeakCountTxt.setText("100");
        maxPeakCountTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                maxPeakCountTxtKeyReleased(evt);
            }
        });

        outputFormatLabel.setText("Output Format");

        outputFormatCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "mzIdentML", "pepXML" }));
        outputFormatCmb.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                outputFormatCmbActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout advancedSettingsPanelLayout = new javax.swing.GroupLayout(advancedSettingsPanel);
        advancedSettingsPanel.setLayout(advancedSettingsPanelLayout);
        advancedSettingsPanelLayout.setHorizontalGroup(
            advancedSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(advancedSettingsPanelLayout.createSequentialGroup()
                .addComponent(useSmartPlusThreeModelLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(useSmartPlusThreeModelCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            .addGroup(advancedSettingsPanelLayout.createSequentialGroup()
                .addComponent(computeXCorrlLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(computeXCorrCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            .addGroup(advancedSettingsPanelLayout.createSequentialGroup()
                .addComponent(ticCutoffPercentageLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(ticCutoffPercentageTxt))
            .addGroup(advancedSettingsPanelLayout.createSequentialGroup()
                .addComponent(numIntensityClassesLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(numIntensityClassesTxt))
            .addGroup(advancedSettingsPanelLayout.createSequentialGroup()
                .addComponent(classSizeMultiplierLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(classSizeMultiplierTxt))
            .addGroup(advancedSettingsPanelLayout.createSequentialGroup()
                .addComponent(numbBatchesLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(numbBatchesTxt))
            .addGroup(advancedSettingsPanelLayout.createSequentialGroup()
                .addComponent(numberMatchesLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(numberMatchesTxt))
            .addGroup(advancedSettingsPanelLayout.createSequentialGroup()
                .addComponent(maxPtmsLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(maxPtmsTxt))
            .addGroup(advancedSettingsPanelLayout.createSequentialGroup()
                .addComponent(numberTerminiLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(terminiCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            .addGroup(advancedSettingsPanelLayout.createSequentialGroup()
                .addComponent(peptideLengthLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(minPepLengthTxt, javax.swing.GroupLayout.DEFAULT_SIZE, 77, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(peptideLengthDividerLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 27, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(maxPepLengthTxt, javax.swing.GroupLayout.DEFAULT_SIZE, 77, Short.MAX_VALUE))
            .addGroup(advancedSettingsPanelLayout.createSequentialGroup()
                .addComponent(precursorMassLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(minPrecursorMassTxt, javax.swing.GroupLayout.DEFAULT_SIZE, 77, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(precursorMassDividerLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 27, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(maxPrecursorMassTxt, javax.swing.GroupLayout.DEFAULT_SIZE, 77, Short.MAX_VALUE))
            .addGroup(advancedSettingsPanelLayout.createSequentialGroup()
                .addComponent(fragmentationMethodLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(fragmentationMethodCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            .addGroup(advancedSettingsPanelLayout.createSequentialGroup()
                .addComponent(maxPeakCountLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(maxPeakCountTxt))
            .addGroup(advancedSettingsPanelLayout.createSequentialGroup()
                .addComponent(outputFormatLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(outputFormatCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        advancedSettingsPanelLayout.setVerticalGroup(
            advancedSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(advancedSettingsPanelLayout.createSequentialGroup()
                .addGroup(advancedSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(minPepLengthTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(maxPepLengthTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(peptideLengthDividerLabel)
                    .addComponent(peptideLengthLabel))
                .addGap(0, 0, 0)
                .addGroup(advancedSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(minPrecursorMassTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(maxPrecursorMassTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(precursorMassDividerLabel)
                    .addComponent(precursorMassLabel))
                .addGap(0, 0, 0)
                .addGroup(advancedSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(numberMatchesLabel)
                    .addComponent(numberMatchesTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(advancedSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(maxPtmsLabel)
                    .addComponent(maxPtmsTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(advancedSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fragmentationMethodCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(fragmentationMethodLabel))
                .addGap(0, 0, 0)
                .addGroup(advancedSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(numberTerminiLabel)
                    .addComponent(terminiCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(advancedSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(useSmartPlusThreeModelLabel)
                    .addComponent(useSmartPlusThreeModelCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(advancedSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(computeXCorrlLabel)
                    .addComponent(computeXCorrCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(advancedSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(ticCutoffPercentageLabel)
                    .addComponent(ticCutoffPercentageTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(advancedSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(numIntensityClassesLabel)
                    .addComponent(numIntensityClassesTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(advancedSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(classSizeMultiplierLabel)
                    .addComponent(classSizeMultiplierTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(advancedSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(numbBatchesLabel)
                    .addComponent(numbBatchesTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(advancedSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(maxPeakCountLabel)
                    .addComponent(maxPeakCountTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, 0)
                .addGroup(advancedSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(outputFormatLabel)
                    .addComponent(outputFormatCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap())
        );

        scrollPane.setViewportView(advancedSettingsPanel);

        javax.swing.GroupLayout advancedSearchSettingsPanelLayout = new javax.swing.GroupLayout(advancedSearchSettingsPanel);
        advancedSearchSettingsPanel.setLayout(advancedSearchSettingsPanelLayout);
        advancedSearchSettingsPanelLayout.setHorizontalGroup(
            advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(advancedSearchSettingsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(scrollPane)
                .addContainerGap())
        );
        advancedSearchSettingsPanelLayout.setVerticalGroup(
            advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(advancedSearchSettingsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(scrollPane)
                .addContainerGap())
        );

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

        advancedSettingsWarningLabel.setText("Click to open the MyriMatch help page.");

        javax.swing.GroupLayout backgroundPanelLayout = new javax.swing.GroupLayout(backgroundPanel);
        backgroundPanel.setLayout(backgroundPanelLayout);
        backgroundPanelLayout.setHorizontalGroup(
            backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, backgroundPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(advancedSearchSettingsPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(backgroundPanelLayout.createSequentialGroup()
                        .addGap(10, 10, 10)
                        .addComponent(openDialogHelpJButton)
                        .addGap(18, 18, 18)
                        .addComponent(advancedSettingsWarningLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(okButton, javax.swing.GroupLayout.PREFERRED_SIZE, 59, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(closeButton)))
                .addContainerGap())
        );
        backgroundPanelLayout.setVerticalGroup(
            backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(backgroundPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(advancedSearchSettingsPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
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
     * Open the MyriMatch help page.
     *
     * @param evt
     */
    private void openDialogHelpJButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_openDialogHelpJButtonActionPerformed
        setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        BareBonesBrowserLaunch.openURL("http://htmlpreview.github.io/?https://github.com/ProteoWizard/pwiz/blob/master/pwiz_tools/Bumbershoot/myrimatch/doc/index.html");
        setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }//GEN-LAST:event_openDialogHelpJButtonActionPerformed

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
    private void numberMatchesTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_numberMatchesTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_numberMatchesTxtKeyReleased

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
    private void minPrecursorMassTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_minPrecursorMassTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_minPrecursorMassTxtKeyReleased

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
    private void ticCutoffPercentageTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_ticCutoffPercentageTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_ticCutoffPercentageTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void numIntensityClassesTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_numIntensityClassesTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_numIntensityClassesTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void classSizeMultiplierTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_classSizeMultiplierTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_classSizeMultiplierTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void numbBatchesTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_numbBatchesTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_numbBatchesTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void maxPeakCountTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_maxPeakCountTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_maxPeakCountTxtKeyReleased

    /**
     * Show warning if pepXML is selected.
     *
     * @param evt
     */
    private void outputFormatCmbActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_outputFormatCmbActionPerformed
        if (outputFormatCmb.getSelectedIndex() != 0 && this.isVisible()) {
            // invoke later to give time for components to update
            SwingUtilities.invokeLater(new Runnable() {
                public void run() {
                    JOptionPane.showMessageDialog(MyriMatchParametersDialog.this, JOptionEditorPane.getJOptionEditorPane(
                            "Note that the MyriMatch pepXML format is not compatible with <a href=\"https://compomics.github.io/projects/peptide-shaker.html\">PeptideShaker</a>."),
                            "Format Warning", JOptionPane.WARNING_MESSAGE);
                }
            });
        }
    }//GEN-LAST:event_outputFormatCmbActionPerformed

    /**
     * Inspects the parameters validity.
     *
     * @param showMessage if true an error messages are shown to the users
     * @return a boolean indicating if the parameters are valid
     */
    public boolean validateInput(boolean showMessage) {

        boolean valid = true;

        valid = GuiUtilities.validateIntegerInput(this, peptideLengthLabel, minPepLengthTxt, "minimum peptide length", "Peptide Length Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, peptideLengthLabel, maxPepLengthTxt, "maximum peptide length", "Peptide Length Error", true, showMessage, valid);
        valid = GuiUtilities.validateDoubleInput(this, precursorMassLabel, minPrecursorMassTxt, "minimum precursor mass", "Precursor Mass Error", true, showMessage, valid);
        valid = GuiUtilities.validateDoubleInput(this, precursorMassLabel, maxPrecursorMassTxt, "maximum precursor mass", "Precursor Mass Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, numberMatchesLabel, numberMatchesTxt, "number of spectrum matches", "Number Spectrum Matches Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, maxPtmsLabel, maxPtmsTxt, "max number of PTMs per peptide", "Peptide PTM Error", true, showMessage, valid);
        valid = GuiUtilities.validateDoubleInput(this, ticCutoffPercentageLabel, ticCutoffPercentageTxt, "TIC cutoff precentage", "TIC Cutoff Percentage Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, numIntensityClassesLabel, numIntensityClassesTxt, "number of intensity classes", "Intensity Classes Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, classSizeMultiplierLabel, classSizeMultiplierTxt, "class size multiplier", "Class Size Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, numbBatchesLabel, numbBatchesTxt, "number of batches", "Number of Batches Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, maxPeakCountLabel, maxPeakCountTxt, "maximum peak count", "Max Peak Count Error", true, showMessage, valid);

        // check that the tic cuttoff is between 0 and 1
        try {
            double temp = Double.parseDouble(ticCutoffPercentageTxt.getText().trim());

            if (temp < 0 || temp > 1) {
                if (showMessage && valid) {
                    JOptionPane.showMessageDialog(this, "The TIC cutoff percentage has to be between 0.0 and 1.0.",
                            "TIC Cutoff Percentage Error", JOptionPane.WARNING_MESSAGE);
                }
                valid = false;
                ticCutoffPercentageLabel.setForeground(Color.RED);
                ticCutoffPercentageLabel.setToolTipText("Please select a valid TIC cutoff percentage [0.0-1.0]");
            }

        } catch (NumberFormatException e) {
            // ignore, handled above
        }

        okButton.setEnabled(valid);

        return valid;
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JPanel advancedSearchSettingsPanel;
    private javax.swing.JPanel advancedSettingsPanel;
    private javax.swing.JLabel advancedSettingsWarningLabel;
    private javax.swing.JPanel backgroundPanel;
    private javax.swing.JLabel classSizeMultiplierLabel;
    private javax.swing.JTextField classSizeMultiplierTxt;
    private javax.swing.JButton closeButton;
    private javax.swing.JComboBox computeXCorrCmb;
    private javax.swing.JLabel computeXCorrlLabel;
    private javax.swing.JComboBox fragmentationMethodCmb;
    private javax.swing.JLabel fragmentationMethodLabel;
    private javax.swing.JLabel maxPeakCountLabel;
    private javax.swing.JTextField maxPeakCountTxt;
    private javax.swing.JTextField maxPepLengthTxt;
    private javax.swing.JTextField maxPrecursorMassTxt;
    private javax.swing.JLabel maxPtmsLabel;
    private javax.swing.JTextField maxPtmsTxt;
    private javax.swing.JTextField minPepLengthTxt;
    private javax.swing.JTextField minPrecursorMassTxt;
    private javax.swing.JLabel numIntensityClassesLabel;
    private javax.swing.JTextField numIntensityClassesTxt;
    private javax.swing.JLabel numbBatchesLabel;
    private javax.swing.JTextField numbBatchesTxt;
    private javax.swing.JLabel numberMatchesLabel;
    private javax.swing.JTextField numberMatchesTxt;
    private javax.swing.JLabel numberTerminiLabel;
    private javax.swing.JButton okButton;
    private javax.swing.JButton openDialogHelpJButton;
    private javax.swing.JComboBox outputFormatCmb;
    private javax.swing.JLabel outputFormatLabel;
    private javax.swing.JLabel peptideLengthDividerLabel;
    private javax.swing.JLabel peptideLengthLabel;
    private javax.swing.JLabel precursorMassDividerLabel;
    private javax.swing.JLabel precursorMassLabel;
    private javax.swing.JScrollPane scrollPane;
    private javax.swing.JComboBox terminiCmb;
    private javax.swing.JLabel ticCutoffPercentageLabel;
    private javax.swing.JTextField ticCutoffPercentageTxt;
    private javax.swing.JComboBox useSmartPlusThreeModelCmb;
    private javax.swing.JLabel useSmartPlusThreeModelLabel;
    // End of variables declaration//GEN-END:variables
}
