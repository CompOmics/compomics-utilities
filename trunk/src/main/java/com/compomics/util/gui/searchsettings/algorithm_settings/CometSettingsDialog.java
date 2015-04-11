package com.compomics.util.gui.searchsettings.algorithm_settings;

import com.compomics.util.examples.BareBonesBrowserLaunch;
import com.compomics.util.experiment.identification.identification_parameters.CometParameters;
import com.compomics.util.gui.GuiUtilities;
import java.awt.Color;
import javax.swing.JOptionPane;
import javax.swing.SwingConstants;

/**
 * Dialog for the Comet specific settings.
 *
 * @author Harald Barsnes
 */
public class CometSettingsDialog extends javax.swing.JDialog {

    /**
     * The Comet parameters class containing the information to display.
     */
    private CometParameters cometParameters;
    /**
     * Boolean indicating whether the used canceled the editing.
     */
    private boolean cancelled = false;

    /**
     * Creates a new CometSettingsDialog.
     *
     * @param parent the parent frame
     * @param cometParameters the Comet parameters
     */
    public CometSettingsDialog(java.awt.Frame parent, CometParameters cometParameters) {
        super(parent, true);
        if (cometParameters != null) {
            this.cometParameters = cometParameters;
        } else {
            this.cometParameters = new CometParameters();
        }
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
        removePrecursorPeakCombo.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        enzymeTypeCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        isotopeCorrectionCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        removePrecursorPeakCombo.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        correlationScoreTypeCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        useSparseMatrixCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        removeMethionineCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
        requireVariablePtmCmb.setRenderer(new com.compomics.util.gui.renderers.AlignedListCellRenderer(SwingConstants.CENTER));
    }

    /**
     * Fills the GUI with the information contained in the Comet settings
     * object.
     */
    private void fillGUI() {

        if (cometParameters.getNumberOfSpectrumMatches() != null) {
            numberMatchesTxt.setText(cometParameters.getNumberOfSpectrumMatches() + "");
        }

        if (cometParameters.getMinPeaks() != null) {
            minPeaksTxt.setText(cometParameters.getMinPeaks() + "");
        }

        if (cometParameters.getMinPeakIntensity() != null) {
            minPeakIntensityTxt.setText(cometParameters.getMinPeakIntensity() + "");
        }

        if (cometParameters.getRemovePrecursor() != null) {
            removePrecursorPeakCombo.setSelectedIndex(cometParameters.getRemovePrecursor());
        }

        if (cometParameters.getRemovePrecursorTolerance() != null) {
            removePrecursorPeakToleranceTxt.setText(cometParameters.getRemovePrecursorTolerance() + "");
        }

        if (cometParameters.getLowerClearMzRange() != null) {
            clearMzRangeLowerTxt.setText(cometParameters.getLowerClearMzRange() + "");
        }
        if (cometParameters.getUpperClearMzRange() != null) {
            clearMzRangeUpperTxt.setText(cometParameters.getUpperClearMzRange() + "");
        }

        if (cometParameters.getEnzymeType() != null) {
            if (cometParameters.getEnzymeType() == 1) {
                enzymeTypeCmb.setSelectedIndex(1);
            } else if (cometParameters.getEnzymeType() == 2) {
                enzymeTypeCmb.setSelectedIndex(0);
            } else if (cometParameters.getEnzymeType() == 8) {
                enzymeTypeCmb.setSelectedIndex(2);
            } else if (cometParameters.getEnzymeType() == 9) {
                enzymeTypeCmb.setSelectedIndex(3);
            }
        }

        if (cometParameters.getIsotopeCorrection() != null) {
            isotopeCorrectionCmb.setSelectedIndex(cometParameters.getIsotopeCorrection());
        }

        if (cometParameters.getMinPrecursorMass() != null) {
            minPrecursorMassTxt.setText(cometParameters.getMinPrecursorMass() + "");
        }
        if (cometParameters.getMaxPrecursorMass() != null) {
            maxPrecursorMassTxt.setText(cometParameters.getMaxPrecursorMass() + "");
        }

        if (cometParameters.getMaxFragmentCharge() != null) {
            maxFragmentChargeTxt.setText(cometParameters.getMaxFragmentCharge() + "");
        }

        if (cometParameters.getRemoveMethionine()) {
            removeMethionineCmb.setSelectedIndex(0);
        } else {
            removeMethionineCmb.setSelectedIndex(1);
        }

        if (cometParameters.getBatchSize() != null) {
            batchSizeTxt.setText(cometParameters.getBatchSize() + "");
        }

        if (cometParameters.getMaxVariableMods() != null) {
            maxPtmsTxt.setText(cometParameters.getMaxVariableMods() + "");
        }

        if (cometParameters.getRequireVariableMods()) {
            requireVariablePtmCmb.setSelectedIndex(0);
        } else {
            requireVariablePtmCmb.setSelectedIndex(1);
        }

        if (cometParameters.getTheoreticalFragmentIonsSumOnly()) {
            correlationScoreTypeCmb.setSelectedIndex(1);
        } else {
            correlationScoreTypeCmb.setSelectedIndex(0);
        }

        if (cometParameters.getFragmentBinOffset() != null) {
            fragmentBinOffsetTxt.setText(cometParameters.getFragmentBinOffset() + "");
        }

        if (cometParameters.getUseSparseMatrix()) {
            useSparseMatrixCmb.setSelectedIndex(0);
        } else {
            useSparseMatrixCmb.setSelectedIndex(1);
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
     * Returns the user selection as Comet parameters object.
     *
     * @return the user selection
     */
    public CometParameters getInput() {

        CometParameters result = new CometParameters();

        String input = minPeaksTxt.getText().trim();
        if (!input.equals("")) {
            result.setMinPeaks(new Integer(input));
        }

        input = minPeakIntensityTxt.getText().trim();
        if (!input.equals("")) {
            result.setMinPeakIntensity(new Double(input));
        }

        result.setRemovePrecursor(removePrecursorPeakCombo.getSelectedIndex());

        input = removePrecursorPeakToleranceTxt.getText().trim();
        if (!input.equals("")) {
            result.setRemovePrecursorTolerance(new Double(input));
        }

        input = clearMzRangeLowerTxt.getText().trim();
        if (!input.equals("")) {
            result.setLowerClearMzRange(new Double(input));
        }
        input = clearMzRangeUpperTxt.getText().trim();
        if (!input.equals("")) {
            result.setUpperClearMzRange(new Double(input));
        }

        int selectedIndex = enzymeTypeCmb.getSelectedIndex();
        if (selectedIndex == 0) {
            result.setEnzymeType(2);
        } else if (selectedIndex == 1) {
            result.setEnzymeType(1);
        } else if (selectedIndex == 2) {
            result.setEnzymeType(8);
        } else if (selectedIndex == 3) {
            result.setEnzymeType(9);
        }

        result.setIsotopeCorrection(isotopeCorrectionCmb.getSelectedIndex());

        input = minPrecursorMassTxt.getText().trim();
        if (!input.equals("")) {
            result.setMinPrecursorMass(new Double(input));
        }
        input = maxPrecursorMassTxt.getText().trim();
        if (!input.equals("")) {
            result.setMaxPrecursorMass(new Double(input));
        }

        input = numberMatchesTxt.getText().trim();
        if (!input.equals("")) {
            result.setNumberOfSpectrumMatches(new Integer(input));
        }

        input = maxFragmentChargeTxt.getText().trim();
        if (!input.equals("")) {
            result.setMaxFragmentCharge(new Integer(input));
        }

        result.setRemoveMethionine(removeMethionineCmb.getSelectedIndex() == 0);

        input = batchSizeTxt.getText().trim();
        if (!input.equals("")) {
            result.setBatchSize(new Integer(input));
        }

        input = maxPtmsTxt.getText().trim();
        if (!input.equals("")) {
            result.setMaxVariableMods(new Integer(input));
        }

        result.setRequireVariableMods(requireVariablePtmCmb.getSelectedIndex() == 0);

        result.setTheoreticalFragmentIonsSumOnly(correlationScoreTypeCmb.getSelectedIndex() == 1);

        input = fragmentBinOffsetTxt.getText().trim();
        if (!input.equals("")) {
            result.setFragmentBinOffset(new Double(input));
        }

        input = fragmentBinOffsetTxt.getText().trim();
        if (!input.equals("")) {
            result.setFragmentBinOffset(new Double(input));
        }

        result.setUseSparseMatrix(useSparseMatrixCmb.getSelectedIndex() == 0);

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
        okButton = new javax.swing.JButton();
        closeButton = new javax.swing.JButton();
        openDialogHelpJButton = new javax.swing.JButton();
        advancedSettingsWarningLabel = new javax.swing.JLabel();
        tabbedPane = new javax.swing.JTabbedPane();
        spectrumProcessingPanel = new javax.swing.JPanel();
        minPeaksLbl = new javax.swing.JLabel();
        minPeaksTxt = new javax.swing.JTextField();
        minPeakIntensityLbl = new javax.swing.JLabel();
        minPeakIntensityTxt = new javax.swing.JTextField();
        removePrecursorPeakLabel = new javax.swing.JLabel();
        removePrecursorPeakCombo = new javax.swing.JComboBox();
        removePrecursorPeakToleranceLbl = new javax.swing.JLabel();
        removePrecursorPeakToleranceTxt = new javax.swing.JTextField();
        clearMzRangeLabel = new javax.swing.JLabel();
        clearMzRangeLowerTxt = new javax.swing.JTextField();
        clearMzRangeDividerLabel = new javax.swing.JLabel();
        clearMzRangeUpperTxt = new javax.swing.JTextField();
        searchSettingsPanel = new javax.swing.JPanel();
        enzymeTypeLabel = new javax.swing.JLabel();
        enzymeTypeCmb = new javax.swing.JComboBox();
        isotopeCorrectionLabel = new javax.swing.JLabel();
        isotopeCorrectionCmb = new javax.swing.JComboBox();
        precursorMassLabel = new javax.swing.JLabel();
        minPrecursorMassTxt = new javax.swing.JTextField();
        precursorMassDividerLabel = new javax.swing.JLabel();
        maxPrecursorMassTxt = new javax.swing.JTextField();
        numberMatchesLabel = new javax.swing.JLabel();
        numberMatchesTxt = new javax.swing.JTextField();
        maxFragmentChargeLabel = new javax.swing.JLabel();
        maxFragmentChargeTxt = new javax.swing.JTextField();
        removeMethionineLabel = new javax.swing.JLabel();
        removeMethionineCmb = new javax.swing.JComboBox();
        batchSizeLabel = new javax.swing.JLabel();
        batchSizeTxt = new javax.swing.JTextField();
        maxPtmsLabel = new javax.swing.JLabel();
        maxPtmsTxt = new javax.swing.JTextField();
        requireVariablePtmLabel = new javax.swing.JLabel();
        requireVariablePtmCmb = new javax.swing.JComboBox();
        fragmentIonsPanel = new javax.swing.JPanel();
        correlationScoreTypeLabel = new javax.swing.JLabel();
        correlationScoreTypeCmb = new javax.swing.JComboBox();
        fragmentBinOffsetLabel = new javax.swing.JLabel();
        fragmentBinOffsetTxt = new javax.swing.JTextField();
        useSparseMatrixLabel = new javax.swing.JLabel();
        useSparseMatrixCmb = new javax.swing.JComboBox();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Advanced Comet Settings");
        setResizable(false);

        backgroundPanel.setBackground(new java.awt.Color(230, 230, 230));

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

        advancedSettingsWarningLabel.setText("Click to open the Comet help page.");

        tabbedPane.setBackground(new java.awt.Color(230, 230, 230));
        tabbedPane.setOpaque(true);

        spectrumProcessingPanel.setBackground(new java.awt.Color(230, 230, 230));
        spectrumProcessingPanel.setPreferredSize(new java.awt.Dimension(518, 143));

        minPeaksLbl.setText("Minimum Number of Peaks");

        minPeaksTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        minPeaksTxt.setText("10");
        minPeaksTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                minPeaksTxtKeyReleased(evt);
            }
        });

        minPeakIntensityLbl.setText("Minimal Peak Intensity");

        minPeakIntensityTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        minPeakIntensityTxt.setText("0.0");
        minPeakIntensityTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                minPeakIntensityTxtKeyReleased(evt);
            }
        });

        removePrecursorPeakLabel.setText("Remove Precursor Peak");

        removePrecursorPeakCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "No", "Yes", "Yes + Charge Reduced " }));

        removePrecursorPeakToleranceLbl.setText("Remove Precursor Peak Tolerance (Da)");

        removePrecursorPeakToleranceTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        removePrecursorPeakToleranceTxt.setText("0.0");
        removePrecursorPeakToleranceTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                removePrecursorPeakToleranceTxtKeyReleased(evt);
            }
        });

        clearMzRangeLabel.setText("Clear m/z Range");

        clearMzRangeLowerTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        clearMzRangeLowerTxt.setText("0.0");
        clearMzRangeLowerTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                clearMzRangeLowerTxtKeyReleased(evt);
            }
        });

        clearMzRangeDividerLabel.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        clearMzRangeDividerLabel.setText("-");

        clearMzRangeUpperTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        clearMzRangeUpperTxt.setText("0.0");
        clearMzRangeUpperTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                clearMzRangeUpperTxtKeyReleased(evt);
            }
        });

        javax.swing.GroupLayout spectrumProcessingPanelLayout = new javax.swing.GroupLayout(spectrumProcessingPanel);
        spectrumProcessingPanel.setLayout(spectrumProcessingPanelLayout);
        spectrumProcessingPanelLayout.setHorizontalGroup(
            spectrumProcessingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(spectrumProcessingPanelLayout.createSequentialGroup()
                .addGap(25, 25, 25)
                .addGroup(spectrumProcessingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(spectrumProcessingPanelLayout.createSequentialGroup()
                        .addComponent(removePrecursorPeakToleranceLbl, javax.swing.GroupLayout.DEFAULT_SIZE, 373, Short.MAX_VALUE)
                        .addGap(120, 120, 120))
                    .addGroup(spectrumProcessingPanelLayout.createSequentialGroup()
                        .addGroup(spectrumProcessingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(spectrumProcessingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                .addComponent(minPeaksLbl, javax.swing.GroupLayout.DEFAULT_SIZE, 250, Short.MAX_VALUE)
                                .addComponent(minPeakIntensityLbl, javax.swing.GroupLayout.DEFAULT_SIZE, 250, Short.MAX_VALUE))
                            .addGroup(spectrumProcessingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                                .addComponent(clearMzRangeLabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, 250, Short.MAX_VALUE)
                                .addComponent(removePrecursorPeakLabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, 250, Short.MAX_VALUE)))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(spectrumProcessingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(spectrumProcessingPanelLayout.createSequentialGroup()
                                .addComponent(clearMzRangeLowerTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 71, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addComponent(clearMzRangeDividerLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addComponent(clearMzRangeUpperTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 77, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addComponent(removePrecursorPeakToleranceTxt)
                            .addComponent(minPeaksTxt)
                            .addComponent(minPeakIntensityTxt, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, 209, Short.MAX_VALUE)
                            .addComponent(removePrecursorPeakCombo, javax.swing.GroupLayout.Alignment.TRAILING, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addGap(30, 30, 30))))
        );

        spectrumProcessingPanelLayout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {clearMzRangeLowerTxt, clearMzRangeUpperTxt});

        spectrumProcessingPanelLayout.setVerticalGroup(
            spectrumProcessingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(spectrumProcessingPanelLayout.createSequentialGroup()
                .addGap(25, 25, 25)
                .addGroup(spectrumProcessingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(minPeaksLbl)
                    .addComponent(minPeaksTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(spectrumProcessingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(minPeakIntensityLbl)
                    .addComponent(minPeakIntensityTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(spectrumProcessingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(removePrecursorPeakLabel)
                    .addComponent(removePrecursorPeakCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(spectrumProcessingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(removePrecursorPeakToleranceLbl)
                    .addComponent(removePrecursorPeakToleranceTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(spectrumProcessingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(clearMzRangeLowerTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(clearMzRangeUpperTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(clearMzRangeDividerLabel)
                    .addComponent(clearMzRangeLabel))
                .addContainerGap(194, Short.MAX_VALUE))
        );

        tabbedPane.addTab("Spectrum", spectrumProcessingPanel);

        searchSettingsPanel.setBackground(new java.awt.Color(230, 230, 230));
        searchSettingsPanel.setPreferredSize(new java.awt.Dimension(518, 143));

        enzymeTypeLabel.setText("Enzyme Type");

        enzymeTypeCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Full-enzyme", "Semi-specific", "Unspecific peptide N-term", "Unspecific peptide C-term" }));

        isotopeCorrectionLabel.setText("Isotope Correction");

        isotopeCorrectionCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "No Correction", "-1, 0, +1, +2, and +3", "-8, -4, 0, +4 and +8" }));

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

        numberMatchesLabel.setText("Number of Spectrum Matches");

        numberMatchesTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        numberMatchesTxt.setText("10");
        numberMatchesTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                numberMatchesTxtKeyReleased(evt);
            }
        });

        maxFragmentChargeLabel.setText("Max Fragment Charge");

        maxFragmentChargeTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        maxFragmentChargeTxt.setText("3");
        maxFragmentChargeTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                maxFragmentChargeTxtKeyReleased(evt);
            }
        });

        removeMethionineLabel.setText("Remove Starting Methionine");

        removeMethionineCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));
        removeMethionineCmb.setSelectedIndex(1);

        batchSizeLabel.setText("Spectrum Batch Size");

        batchSizeTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        batchSizeTxt.setText("0");
        batchSizeTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                batchSizeTxtKeyReleased(evt);
            }
        });

        maxPtmsLabel.setText("Max Variable PTMs per Peptide");

        maxPtmsTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        maxPtmsTxt.setText("10");
        maxPtmsTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                maxPtmsTxtKeyReleased(evt);
            }
        });

        requireVariablePtmLabel.setText("Require Variable PTM");

        requireVariablePtmCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));
        requireVariablePtmCmb.setSelectedIndex(1);

        javax.swing.GroupLayout searchSettingsPanelLayout = new javax.swing.GroupLayout(searchSettingsPanel);
        searchSettingsPanel.setLayout(searchSettingsPanelLayout);
        searchSettingsPanelLayout.setHorizontalGroup(
            searchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(searchSettingsPanelLayout.createSequentialGroup()
                .addGap(25, 25, 25)
                .addGroup(searchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(searchSettingsPanelLayout.createSequentialGroup()
                        .addGroup(searchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(searchSettingsPanelLayout.createSequentialGroup()
                                .addGroup(searchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                    .addComponent(precursorMassLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .addGroup(searchSettingsPanelLayout.createSequentialGroup()
                                        .addComponent(maxFragmentChargeLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 250, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addGap(4, 4, 4)))
                                .addGap(0, 0, Short.MAX_VALUE))
                            .addGroup(searchSettingsPanelLayout.createSequentialGroup()
                                .addGroup(searchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addGroup(searchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                        .addComponent(removeMethionineLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(batchSizeLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 250, javax.swing.GroupLayout.PREFERRED_SIZE))
                                    .addComponent(numberMatchesLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 250, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(isotopeCorrectionLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 250, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(enzymeTypeLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 250, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 23, Short.MAX_VALUE)))
                        .addGroup(searchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(numberMatchesTxt)
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, searchSettingsPanelLayout.createSequentialGroup()
                                .addComponent(minPrecursorMassTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 77, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(precursorMassDividerLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 35, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addComponent(maxPrecursorMassTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 83, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addComponent(enzymeTypeCmb, javax.swing.GroupLayout.Alignment.TRAILING, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(maxFragmentChargeTxt)
                            .addComponent(removeMethionineCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(batchSizeTxt, javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(isotopeCorrectionCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                    .addGroup(searchSettingsPanelLayout.createSequentialGroup()
                        .addComponent(maxPtmsLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 250, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(maxPtmsTxt))
                    .addGroup(searchSettingsPanelLayout.createSequentialGroup()
                        .addComponent(requireVariablePtmLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 250, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(requireVariablePtmCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                .addGap(30, 30, 30))
        );

        searchSettingsPanelLayout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {maxPtmsLabel, requireVariablePtmLabel});

        searchSettingsPanelLayout.setVerticalGroup(
            searchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(searchSettingsPanelLayout.createSequentialGroup()
                .addGap(25, 25, 25)
                .addGroup(searchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(enzymeTypeLabel)
                    .addComponent(enzymeTypeCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(searchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(isotopeCorrectionLabel)
                    .addComponent(isotopeCorrectionCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(searchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(minPrecursorMassTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(maxPrecursorMassTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(precursorMassDividerLabel)
                    .addComponent(precursorMassLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(searchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(numberMatchesLabel)
                    .addComponent(numberMatchesTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(searchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(maxFragmentChargeLabel)
                    .addComponent(maxFragmentChargeTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(searchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(removeMethionineLabel)
                    .addComponent(removeMethionineCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(searchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(batchSizeLabel)
                    .addComponent(batchSizeTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(searchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(maxPtmsLabel)
                    .addComponent(maxPtmsTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(searchSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(requireVariablePtmCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(requireVariablePtmLabel))
                .addContainerGap(90, Short.MAX_VALUE))
        );

        tabbedPane.addTab("Search", searchSettingsPanel);

        fragmentIonsPanel.setBackground(new java.awt.Color(230, 230, 230));

        correlationScoreTypeLabel.setText("Correlation Score Type");

        correlationScoreTypeCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Summed Intensities +  Flanking", "Summed Intensities" }));

        fragmentBinOffsetLabel.setText("Fragment Bin Offset");

        fragmentBinOffsetTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        fragmentBinOffsetTxt.setText("0.4");
        fragmentBinOffsetTxt.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                fragmentBinOffsetTxtKeyReleased(evt);
            }
        });

        useSparseMatrixLabel.setText("Use Sparse Matrix");

        useSparseMatrixCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));
        useSparseMatrixCmb.setSelectedIndex(1);

        javax.swing.GroupLayout fragmentIonsPanelLayout = new javax.swing.GroupLayout(fragmentIonsPanel);
        fragmentIonsPanel.setLayout(fragmentIonsPanelLayout);
        fragmentIonsPanelLayout.setHorizontalGroup(
            fragmentIonsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, fragmentIonsPanelLayout.createSequentialGroup()
                .addGap(25, 25, 25)
                .addGroup(fragmentIonsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(correlationScoreTypeLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 250, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(fragmentBinOffsetLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 250, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(useSparseMatrixLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 250, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(fragmentIonsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(useSparseMatrixCmb, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(fragmentBinOffsetTxt)
                    .addComponent(correlationScoreTypeCmb, 0, 209, Short.MAX_VALUE))
                .addGap(30, 30, 30))
        );
        fragmentIonsPanelLayout.setVerticalGroup(
            fragmentIonsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(fragmentIonsPanelLayout.createSequentialGroup()
                .addGap(25, 25, 25)
                .addGroup(fragmentIonsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(correlationScoreTypeLabel)
                    .addComponent(correlationScoreTypeCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(fragmentIonsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fragmentBinOffsetLabel)
                    .addComponent(fragmentBinOffsetTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(fragmentIonsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(useSparseMatrixLabel)
                    .addComponent(useSparseMatrixCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(246, Short.MAX_VALUE))
        );

        tabbedPane.addTab("Fragment Ions", fragmentIonsPanel);

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
                .addComponent(tabbedPane, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        backgroundPanelLayout.setVerticalGroup(
            backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(backgroundPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(tabbedPane, javax.swing.GroupLayout.PREFERRED_SIZE, 371, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
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
     * Open the Comet help page.
     *
     * @param evt
     */
    private void openDialogHelpJButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_openDialogHelpJButtonActionPerformed
        setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        BareBonesBrowserLaunch.openURL("http://comet-ms.sourceforge.net");
        setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }//GEN-LAST:event_openDialogHelpJButtonActionPerformed

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void clearMzRangeLowerTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_clearMzRangeLowerTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_clearMzRangeLowerTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void clearMzRangeUpperTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_clearMzRangeUpperTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_clearMzRangeUpperTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void fragmentBinOffsetTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_fragmentBinOffsetTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_fragmentBinOffsetTxtKeyReleased

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
    private void minPeakIntensityTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_minPeakIntensityTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_minPeakIntensityTxtKeyReleased

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
    private void maxPtmsTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_maxPtmsTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_maxPtmsTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void batchSizeTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_batchSizeTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_batchSizeTxtKeyReleased

    /**
     * Validate the input.
     *
     * @param evt
     */
    private void maxFragmentChargeTxtKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_maxFragmentChargeTxtKeyReleased
        validateInput(false);
    }//GEN-LAST:event_maxFragmentChargeTxtKeyReleased

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
     * Inspects the parameter validity.
     *
     * @param showMessage if true an error messages are shown to the users
     * @return a boolean indicating if the parameters are valid
     */
    public boolean validateInput(boolean showMessage) {

        boolean valid = true;

        valid = GuiUtilities.validateIntegerInput(this, minPeaksLbl, minPeaksTxt, "minimum number of peaks", "Minimum Number of Peaks Error", true, showMessage, valid);
        valid = GuiUtilities.validateDoubleInput(this, minPeakIntensityLbl, minPeakIntensityTxt, "minimum peak intensity", "Minimim Peak Intensity Error", true, showMessage, valid);
        valid = GuiUtilities.validateDoubleInput(this, removePrecursorPeakLabel, removePrecursorPeakToleranceTxt, "remove precursor peak tolerance", "Precursor Peak Error", true, showMessage, valid);
        valid = GuiUtilities.validateDoubleInput(this, clearMzRangeLabel, clearMzRangeLowerTxt, "lower clear mz range", "Clear MZ Range Error", true, showMessage, valid);
        valid = GuiUtilities.validateDoubleInput(this, clearMzRangeLabel, clearMzRangeUpperTxt, "upper clear mz range", "Clear MZ Range Error", true, showMessage, valid);
        valid = GuiUtilities.validateDoubleInput(this, precursorMassLabel, minPrecursorMassTxt, "minimum precursor mass", "Precursor Mass Error", true, showMessage, valid);
        valid = GuiUtilities.validateDoubleInput(this, precursorMassLabel, maxPrecursorMassTxt, "maximum precursor mass", "Precursor Mass Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, numberMatchesLabel, numberMatchesTxt, "number of spectrum matches", "Number of Spectrum Matches Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, maxFragmentChargeLabel, maxFragmentChargeTxt, "maximum fragment charge", "Maximum Fragment Charge Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, batchSizeLabel, batchSizeTxt, "batch size", "Batch Size Error", true, showMessage, valid);
        valid = GuiUtilities.validateIntegerInput(this, maxPtmsLabel, maxPtmsTxt, "maximum number of variable PTMs", "Variable PTMs Error", true, showMessage, valid);
        valid = GuiUtilities.validateDoubleInput(this, fragmentBinOffsetLabel, fragmentBinOffsetTxt, "fragment bin offset", "Fragment Bin Offset Error", true, showMessage, valid);

        // clear range: the low value should be lower than the high value
        try {
            double lowValue = Double.parseDouble(clearMzRangeLowerTxt.getText().trim());
            double highValue = Double.parseDouble(clearMzRangeUpperTxt.getText().trim());

            if (lowValue > highValue) {
                if (showMessage && valid) {
                    JOptionPane.showMessageDialog(this, "The lower range value has to be smaller than the upper range value.",
                            "Clear MZ Range Error", JOptionPane.WARNING_MESSAGE);
                }
                valid = false;
                clearMzRangeLabel.setForeground(Color.RED);
                clearMzRangeLabel.setToolTipText("Please select a valid range (upper <= higher)");
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

        // @TODO: add more tests?
        okButton.setEnabled(valid);

        return valid;
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel advancedSettingsWarningLabel;
    private javax.swing.JPanel backgroundPanel;
    private javax.swing.JLabel batchSizeLabel;
    private javax.swing.JTextField batchSizeTxt;
    private javax.swing.JLabel clearMzRangeDividerLabel;
    private javax.swing.JLabel clearMzRangeLabel;
    private javax.swing.JTextField clearMzRangeLowerTxt;
    private javax.swing.JTextField clearMzRangeUpperTxt;
    private javax.swing.JButton closeButton;
    private javax.swing.JComboBox correlationScoreTypeCmb;
    private javax.swing.JLabel correlationScoreTypeLabel;
    private javax.swing.JComboBox enzymeTypeCmb;
    private javax.swing.JLabel enzymeTypeLabel;
    private javax.swing.JLabel fragmentBinOffsetLabel;
    private javax.swing.JTextField fragmentBinOffsetTxt;
    private javax.swing.JPanel fragmentIonsPanel;
    private javax.swing.JComboBox isotopeCorrectionCmb;
    private javax.swing.JLabel isotopeCorrectionLabel;
    private javax.swing.JLabel maxFragmentChargeLabel;
    private javax.swing.JTextField maxFragmentChargeTxt;
    private javax.swing.JTextField maxPrecursorMassTxt;
    private javax.swing.JLabel maxPtmsLabel;
    private javax.swing.JTextField maxPtmsTxt;
    private javax.swing.JLabel minPeakIntensityLbl;
    private javax.swing.JTextField minPeakIntensityTxt;
    private javax.swing.JLabel minPeaksLbl;
    private javax.swing.JTextField minPeaksTxt;
    private javax.swing.JTextField minPrecursorMassTxt;
    private javax.swing.JLabel numberMatchesLabel;
    private javax.swing.JTextField numberMatchesTxt;
    private javax.swing.JButton okButton;
    private javax.swing.JButton openDialogHelpJButton;
    private javax.swing.JLabel precursorMassDividerLabel;
    private javax.swing.JLabel precursorMassLabel;
    private javax.swing.JComboBox removeMethionineCmb;
    private javax.swing.JLabel removeMethionineLabel;
    private javax.swing.JComboBox removePrecursorPeakCombo;
    private javax.swing.JLabel removePrecursorPeakLabel;
    private javax.swing.JLabel removePrecursorPeakToleranceLbl;
    private javax.swing.JTextField removePrecursorPeakToleranceTxt;
    private javax.swing.JComboBox requireVariablePtmCmb;
    private javax.swing.JLabel requireVariablePtmLabel;
    private javax.swing.JPanel searchSettingsPanel;
    private javax.swing.JPanel spectrumProcessingPanel;
    private javax.swing.JTabbedPane tabbedPane;
    private javax.swing.JComboBox useSparseMatrixCmb;
    private javax.swing.JLabel useSparseMatrixLabel;
    // End of variables declaration//GEN-END:variables
}
