package com.compomics.util.gui.searchsettings.algorithm_settings;

import com.compomics.util.examples.BareBonesBrowserLaunch;
import com.compomics.util.experiment.identification.identification_parameters.MsgfParameters;
import com.compomics.util.gui.GuiUtilities;
import java.awt.Color;
import javax.swing.JOptionPane;
import javax.swing.SwingConstants;

/**
 * Dialog for the MS-GF+ specific settings.
 *
 * @author Harald Barsnes
 */
public class MsgfSettingsDialog extends javax.swing.JDialog {

    /**
     * The MS-GF+ parameters class containing the information to display.
     */
    private MsgfParameters msgfParameters;
    /**
     * Boolean indicating whether the used canceled the editing.
     */
    private boolean cancelled = false;

    /**
     * Creates new form MsgfSettingsDialog.
     *
     * @param parent the parent frame
     * @param msgfParameters the MS-GF+ parameters
     */
    public MsgfSettingsDialog(java.awt.Frame parent, MsgfParameters msgfParameters) {
        super(parent, true);
        this.msgfParameters = msgfParameters;
        initComponents();
        setUpGUI();
        fillGUI();
        validateInput(false);
        setLocationRelativeTo(parent);
        setVisible(true);
    }

    /**
     * Set up the GUI.
     */
    private void setUpGUI() {
        decoyDatabaseCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        instrumentCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        fragmentationMethodCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        protocolCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        additionalOutputCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
    }

    /**
     * Fills the GUI with the information contained in the MS-GF+ settings
     * object.
     */
    private void fillGUI() {

        if (msgfParameters.searchDecoyDatabase()) {
            decoyDatabaseCmb.setSelectedIndex(0);
        } else {
            decoyDatabaseCmb.setSelectedIndex(1);
        }

        instrumentCmb.setSelectedIndex(msgfParameters.getInstrumentID());
        fragmentationMethodCmb.setSelectedIndex(msgfParameters.getFragmentationType());
        protocolCmb.setSelectedIndex(msgfParameters.getProtocol());

        if (msgfParameters.getMinPeptideLength() != null) {
            minPepLengthTxt.setText(msgfParameters.getMinPeptideLength() + "");
        }
        if (msgfParameters.getMaxPeptideLength() != null) {
            maxPepLengthTxt.setText(msgfParameters.getMaxPeptideLength() + "");
        }
        if (msgfParameters.getNumberOfSpectrumMarches() != null) {
            numberMatchesTxt.setText(msgfParameters.getNumberOfSpectrumMarches() + "");
        }

        if (msgfParameters.isAdditionalOutput()) {
            additionalOutputCmb.setSelectedIndex(0);
        } else {
            additionalOutputCmb.setSelectedIndex(1);
        }

        if (msgfParameters.getLowerIsotopeErrorRange() != null) {
            lowIsotopeErrorRangeTxt.setText(msgfParameters.getLowerIsotopeErrorRange() + "");
        }
        if (msgfParameters.getUpperIsotopeErrorRange() != null) {
            highIsotopeErrorRangeTxt.setText(msgfParameters.getUpperIsotopeErrorRange() + "");
        }
        if (msgfParameters.getMaxEValue() != null) {
            eValueCutoffTxt.setText(msgfParameters.getMaxEValue() + "");
        }
        if (msgfParameters.getMaxEValue() != null) {
            eValueCutoffTxt.setText(msgfParameters.getMaxEValue() + "");
        }
        if (msgfParameters.getNumberTolerableTermini() != null) {
            numberTerminiTxt.setText(msgfParameters.getNumberTolerableTermini() + "");
        }
        if (msgfParameters.getNumberOfPtmsPerPeptide() != null) {
            maxPtmsTxt.setText(msgfParameters.getNumberOfPtmsPerPeptide() + "");
        }
    }

    /**
     * Indicates whether the user canceled the process.
     *
     * @return true if cancel was pressed
     */
    public boolean isCancelled() {
        return cancelled;
    }

    /**
     * Returns the user selection as MS-GF+ parameters object.
     *
     * @return the user selection
     */
    public MsgfParameters getInput() {

        MsgfParameters result = new MsgfParameters();

        result.setSearchDecoyDatabase(decoyDatabaseCmb.getSelectedIndex() == 0);
        result.setInstrumentID(instrumentCmb.getSelectedIndex());
        result.setFragmentationType(fragmentationMethodCmb.getSelectedIndex());
        result.setProtocol(protocolCmb.getSelectedIndex());

        String input = minPepLengthTxt.getText().trim();
        if (!input.equals("")) {
            result.setMinPeptideLength(new Integer(input));
        }
        input = maxPepLengthTxt.getText().trim();
        if (!input.equals("")) {
            result.setMaxPeptideLength(new Integer(input));
        }
        input = numberMatchesTxt.getText().trim();
        if (!input.equals("")) {
            result.setNumberOfSpectrumMarches(new Integer(input));
        }

        result.setAdditionalOutput(additionalOutputCmb.getSelectedIndex() == 0);

        input = lowIsotopeErrorRangeTxt.getText().trim();
        if (!input.equals("")) {
            result.setLowerIsotopeErrorRange(new Integer(input));
        }
        input = highIsotopeErrorRangeTxt.getText().trim();
        if (!input.equals("")) {
            result.setUpperIsotopeErrorRange(new Integer(input));
        }
        input = eValueCutoffTxt.getText().trim();
        if (!input.equals("")) {
            result.setMaxEValue(new Double(input));
        }
        input = numberTerminiTxt.getText().trim();
        if (!input.equals("")) {
            result.setNumberTolerableTermini(new Integer(input));
        }
        input = maxPtmsTxt.getText().trim();
        if (!input.equals("")) {
            result.setNumberOfPtmsPerPeptide(new Integer(input));
        }

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
        instrumentCmb = new javax.swing.JComboBox();
        decoyDatabaseCmb = new javax.swing.JComboBox();
        decoyDatabaseLabel = new javax.swing.JLabel();
        instrumentLabel = new javax.swing.JLabel();
        fragmentationMethodLabel = new javax.swing.JLabel();
        fragmentationMethodCmb = new javax.swing.JComboBox();
        protocolLabel = new javax.swing.JLabel();
        protocolCmb = new javax.swing.JComboBox();
        minPepLengthTxt = new javax.swing.JTextField();
        peptideLengthDividerLabel = new javax.swing.JLabel();
        maxPepLengthTxt = new javax.swing.JTextField();
        peptideLengthLabel = new javax.swing.JLabel();
        numberMatchesLabel = new javax.swing.JLabel();
        numberMatchesTxt = new javax.swing.JTextField();
        additionalOutputlLabel = new javax.swing.JLabel();
        additionalOutputCmb = new javax.swing.JComboBox();
        isotopeErrorRangeLabel = new javax.swing.JLabel();
        lowIsotopeErrorRangeTxt = new javax.swing.JTextField();
        highIsotopeErrorRangeTxt = new javax.swing.JTextField();
        isotopeErrorRangeDividerLabel = new javax.swing.JLabel();
        eValueCutoffLabel = new javax.swing.JLabel();
        eValueCutoffTxt = new javax.swing.JTextField();
        numberTerminiLabel = new javax.swing.JLabel();
        numberTerminiTxt = new javax.swing.JTextField();
        maxPtmsLabel = new javax.swing.JLabel();
        maxPtmsTxt = new javax.swing.JTextField();
        okButton = new javax.swing.JButton();
        closeButton = new javax.swing.JButton();
        openDialogHelpJButton = new javax.swing.JButton();
        advancedSettingsWarningLabel = new javax.swing.JLabel();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Advanced MS-GF+ Settings");
        setResizable(false);

        backgroundPanel.setBackground(new java.awt.Color(230, 230, 230));

        advancedSearchSettingsPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Search Settings"));
        advancedSearchSettingsPanel.setOpaque(false);

        instrumentCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Low-res LCQ/LTQ", "Orbitrap/FTICR", "TOF", "Q-Exactive" }));

        decoyDatabaseCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));
        decoyDatabaseCmb.setSelectedIndex(1);

        decoyDatabaseLabel.setText("Search Decoy Database");

        instrumentLabel.setText("MS/MS Detector");

        fragmentationMethodLabel.setText("Fragmentation Method");

        fragmentationMethodCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Automatic", "CID", "ETD", "HCD" }));

        protocolLabel.setText("Protocol");

        protocolCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "No Protocol", "Phosphorylation", "iTRAQ", "iTRAQPhospho", "TMT" }));

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
        maxPepLengthTxt.setText("40");
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

        additionalOutputlLabel.setText("Additional Output");

        additionalOutputCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));
        additionalOutputCmb.setSelectedIndex(1);

        isotopeErrorRangeLabel.setText("Isotope Error Range");

        lowIsotopeErrorRangeTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        lowIsotopeErrorRangeTxt.setText("0");
        lowIsotopeErrorRangeTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                lowIsotopeErrorRangeTxtKeyReleased(evt);
            }
        });

        highIsotopeErrorRangeTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        highIsotopeErrorRangeTxt.setText("1");
        highIsotopeErrorRangeTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                highIsotopeErrorRangeTxtKeyReleased(evt);
            }
        });

        isotopeErrorRangeDividerLabel.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        isotopeErrorRangeDividerLabel.setText("-");

        eValueCutoffLabel.setText("E-value Cutoff");

        eValueCutoffTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        eValueCutoffTxt.setText("100");
        eValueCutoffTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                eValueCutoffTxtKeyReleased(evt);
            }
        });

        numberTerminiLabel.setText("Number of Tolerable Termini");

        numberTerminiTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        numberTerminiTxt.setText("2");
        numberTerminiTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                numberTerminiTxtKeyReleased(evt);
            }
        });

        maxPtmsLabel.setText("Max Variable PTMs per Peptide");

        maxPtmsTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        maxPtmsTxt.setText("2");
        maxPtmsTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                maxPtmsTxtKeyReleased(evt);
            }
        });

        javax.swing.GroupLayout advancedSearchSettingsPanelLayout = new javax.swing.GroupLayout(advancedSearchSettingsPanel);
        advancedSearchSettingsPanel.setLayout(advancedSearchSettingsPanelLayout);
        advancedSearchSettingsPanelLayout.setHorizontalGroup(
            advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(advancedSearchSettingsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(advancedSearchSettingsPanelLayout.createSequentialGroup()
                        .addComponent(numberMatchesLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(numberMatchesTxt))
                    .addGroup(advancedSearchSettingsPanelLayout.createSequentialGroup()
                        .addGroup(advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(instrumentLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(decoyDatabaseLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(fragmentationMethodLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(fragmentationMethodCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(instrumentCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(decoyDatabaseCmb, javax.swing.GroupLayout.Alignment.TRAILING, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                    .addGroup(advancedSearchSettingsPanelLayout.createSequentialGroup()
                        .addGroup(advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(protocolLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(peptideLengthLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(advancedSearchSettingsPanelLayout.createSequentialGroup()
                                .addComponent(minPepLengthTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 71, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(peptideLengthDividerLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 27, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(maxPepLengthTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 77, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addComponent(protocolCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                    .addGroup(advancedSearchSettingsPanelLayout.createSequentialGroup()
                        .addComponent(additionalOutputlLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(additionalOutputCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addGroup(advancedSearchSettingsPanelLayout.createSequentialGroup()
                        .addComponent(isotopeErrorRangeLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(lowIsotopeErrorRangeTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 76, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(isotopeErrorRangeDividerLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 11, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(highIsotopeErrorRangeTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 76, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(advancedSearchSettingsPanelLayout.createSequentialGroup()
                        .addComponent(eValueCutoffLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(eValueCutoffTxt))
                    .addGroup(advancedSearchSettingsPanelLayout.createSequentialGroup()
                        .addComponent(numberTerminiLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(numberTerminiTxt))
                    .addGroup(advancedSearchSettingsPanelLayout.createSequentialGroup()
                        .addComponent(maxPtmsLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 200, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(maxPtmsTxt)))
                .addContainerGap())
        );

        advancedSearchSettingsPanelLayout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {highIsotopeErrorRangeTxt, lowIsotopeErrorRangeTxt, maxPepLengthTxt, minPepLengthTxt});

        advancedSearchSettingsPanelLayout.setVerticalGroup(
            advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(advancedSearchSettingsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(decoyDatabaseLabel)
                    .addComponent(decoyDatabaseCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(instrumentCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(instrumentLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fragmentationMethodCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(fragmentationMethodLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(protocolCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(protocolLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(minPepLengthTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(maxPepLengthTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(peptideLengthDividerLabel)
                    .addComponent(peptideLengthLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(numberMatchesLabel)
                    .addComponent(numberMatchesTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(additionalOutputlLabel)
                    .addComponent(additionalOutputCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(lowIsotopeErrorRangeTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(highIsotopeErrorRangeTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(isotopeErrorRangeDividerLabel)
                    .addComponent(isotopeErrorRangeLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(numberTerminiLabel)
                    .addComponent(numberTerminiTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(maxPtmsLabel)
                    .addComponent(maxPtmsTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(advancedSearchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(eValueCutoffLabel)
                    .addComponent(eValueCutoffTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
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

        advancedSettingsWarningLabel.setText("Click to open the MS-GF+ help page.");

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
                .addComponent(advancedSearchSettingsPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
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
            .addComponent(backgroundPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
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
     * Open the MS-GF+ help page.
     *
     * @param evt
     */
    private void openDialogHelpJButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_openDialogHelpJButtonActionPerformed
        setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        BareBonesBrowserLaunch.openURL("https://bix-lab.ucsd.edu/pages/viewpage.action?pageId=13533355");
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
    private void lowIsotopeErrorRangeTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_lowIsotopeErrorRangeTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_lowIsotopeErrorRangeTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void highIsotopeErrorRangeTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_highIsotopeErrorRangeTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_highIsotopeErrorRangeTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void eValueCutoffTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_eValueCutoffTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_eValueCutoffTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void numberTerminiTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_numberTerminiTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_numberTerminiTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void maxPtmsTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_maxPtmsTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_maxPtmsTxtKeyReleased

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
        valid = GuiUtilities.validateIntegerInput(this, numberMatchesLabel, numberMatchesTxt, "number of spectrum matches", "Number Spectrum Matches Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, isotopeErrorRangeLabel, lowIsotopeErrorRangeTxt, "lower isotope range", "Isotope Range Error", false, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, isotopeErrorRangeLabel, highIsotopeErrorRangeTxt, "upper isotope range", "Isotope Range Error", false, showMessage, valid);
        valid = GuiUtilities.validateDoubleInput(this, eValueCutoffLabel, eValueCutoffTxt, "e-value cuttoff limit", "E-value Cutoff Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, numberTerminiLabel, numberTerminiTxt, "number of tolerable termini", "Peptide Termini Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, maxPtmsLabel, maxPtmsTxt, "max number of PTMs per peptide", "Peptide PTM Error", true, showMessage, valid);

        // isotope range: the low value should be lower than the high value
        try {
            int lowValue = Integer.parseInt(lowIsotopeErrorRangeTxt.getText().trim());
            int highValue = Integer.parseInt(highIsotopeErrorRangeTxt.getText().trim());

            if (lowValue > highValue) {
                if (showMessage && valid) {
                    JOptionPane.showMessageDialog(this, "The lower range value has to be smaller than the upper range value.",
                            "Isotope Range Error", JOptionPane.WARNING_MESSAGE);
                }
                valid = false;
                isotopeErrorRangeLabel.setForeground(Color.RED);
                isotopeErrorRangeLabel.setToolTipText("Please select a valid range (upper < higher))");
            }
        } catch (NumberFormatException e) {
            // ignore, handled above
        }
        
        okButton.setEnabled(valid);

        return valid;
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JComboBox additionalOutputCmb;
    private javax.swing.JLabel additionalOutputlLabel;
    private javax.swing.JPanel advancedSearchSettingsPanel;
    private javax.swing.JLabel advancedSettingsWarningLabel;
    private javax.swing.JPanel backgroundPanel;
    private javax.swing.JButton closeButton;
    private javax.swing.JComboBox decoyDatabaseCmb;
    private javax.swing.JLabel decoyDatabaseLabel;
    private javax.swing.JLabel eValueCutoffLabel;
    private javax.swing.JTextField eValueCutoffTxt;
    private javax.swing.JComboBox fragmentationMethodCmb;
    private javax.swing.JLabel fragmentationMethodLabel;
    private javax.swing.JTextField highIsotopeErrorRangeTxt;
    private javax.swing.JComboBox instrumentCmb;
    private javax.swing.JLabel instrumentLabel;
    private javax.swing.JLabel isotopeErrorRangeDividerLabel;
    private javax.swing.JLabel isotopeErrorRangeLabel;
    private javax.swing.JTextField lowIsotopeErrorRangeTxt;
    private javax.swing.JTextField maxPepLengthTxt;
    private javax.swing.JLabel maxPtmsLabel;
    private javax.swing.JTextField maxPtmsTxt;
    private javax.swing.JTextField minPepLengthTxt;
    private javax.swing.JLabel numberMatchesLabel;
    private javax.swing.JTextField numberMatchesTxt;
    private javax.swing.JLabel numberTerminiLabel;
    private javax.swing.JTextField numberTerminiTxt;
    private javax.swing.JButton okButton;
    private javax.swing.JButton openDialogHelpJButton;
    private javax.swing.JLabel peptideLengthDividerLabel;
    private javax.swing.JLabel peptideLengthLabel;
    private javax.swing.JComboBox protocolCmb;
    private javax.swing.JLabel protocolLabel;
    // End of variables declaration//GEN-END:variables
}
