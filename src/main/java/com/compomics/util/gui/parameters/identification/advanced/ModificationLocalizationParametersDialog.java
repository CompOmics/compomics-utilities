package com.compomics.util.gui.parameters.identification.advanced;

import com.compomics.util.experiment.identification.modification.ModificationLocalizationScore;
import com.compomics.util.gui.error_handlers.HelpDialog;
import com.compomics.util.gui.renderers.AlignedListCellRenderer;
import com.compomics.util.parameters.identification.advanced.ModificationLocalizationParameters;
import java.awt.Dialog;
import java.awt.Toolkit;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JOptionPane;
import javax.swing.SwingConstants;

/**
 * PTMLocalizationParametersDialog.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class ModificationLocalizationParametersDialog extends javax.swing.JDialog {

    /**
     * Empty default constructor
     */
    public ModificationLocalizationParametersDialog() {
    }

    /**
     * The parent frame.
     */
    private java.awt.Frame parentFrame;
    /**
     * Boolean indicating whether the user canceled the editing.
     */
    private boolean canceled = false;
    /**
     * Boolean indicating whether the settings can be edited by the user.
     */
    private boolean editable;

    /**
     * Creates a new PTMLocalizationParametersDialog with a frame as owner.
     *
     * @param parentFrame a parent frame
     * @param ptmScoringPreferences the PTM scoring preferences to display
     * @param editable boolean indicating whether the settings can be edited by
     * the user
     */
    public ModificationLocalizationParametersDialog(java.awt.Frame parentFrame, ModificationLocalizationParameters ptmScoringPreferences, boolean editable) {
        super(parentFrame, true);
        this.parentFrame = parentFrame;
        this.editable = editable;
        initComponents();
        setUpGui();
        populateGUI(ptmScoringPreferences);
        setLocationRelativeTo(parentFrame);
        setVisible(true);
    }

    /**
     * Creates a new PTMLocalizationParametersDialog with a dialog as owner.
     *
     * @param owner the dialog owner
     * @param parentFrame a parent frame
     * @param ptmScoringPreferences the PTM scoring preferences to display
     * @param editable boolean indicating whether the settings can be edited by
     * the user
     */
    public ModificationLocalizationParametersDialog(Dialog owner, java.awt.Frame parentFrame, ModificationLocalizationParameters ptmScoringPreferences, boolean editable) {
        super(owner, true);
        this.parentFrame = parentFrame;
        this.editable = editable;
        initComponents();
        setUpGui();
        populateGUI(ptmScoringPreferences);
        setLocationRelativeTo(owner);
        setVisible(true);
    }

    /**
     * Set up the GUI.
     */
    private void setUpGui() {

        neutralLossesCmb.setRenderer(new AlignedListCellRenderer(SwingConstants.CENTER));
        alignOnConfidentCmb.setRenderer(new AlignedListCellRenderer(SwingConstants.CENTER));
        scoreCmb.setRenderer(new AlignedListCellRenderer(SwingConstants.CENTER));

        scoreCmb.setEnabled(editable);
        neutralLossesCmb.setEnabled(editable);
        alignOnConfidentCmb.setEnabled(editable);
        thresholdTxt.setEditable(editable);
        thresholdTxt.setEnabled(editable);
    }

    /**
     * Fills the GUI with the given settings.
     *
     * @param spectrumCountingPreferences the PTM scoring preferences to display
     */
    private void populateGUI(ModificationLocalizationParameters ptmScoringPreferences) {

        scoreCmb.setSelectedItem(ptmScoringPreferences.getSelectedProbabilisticScore());

        if (ptmScoringPreferences.isProbabilisticScoreCalculation()) {
            scoreCmb.setEnabled(editable);
            neutralLossesCmb.setEnabled(editable);
            if (ptmScoringPreferences.isProbabilisticScoreNeutralLosses()) {
                neutralLossesCmb.setSelectedIndex(0);
            } else {
                neutralLossesCmb.setSelectedIndex(1);
            }

            thresholdTxt.setEnabled(editable);
            thresholdTxt.setEditable(editable);
            thresholdTxt.setText(ptmScoringPreferences.getProbabilisticScoreThreshold() + "");

        } else {

            neutralLossesCmb.setEnabled(false);
            thresholdTxt.setEnabled(false);
            thresholdTxt.setEditable(false);

        }

        if (ptmScoringPreferences.getAlignNonConfidentModifications()) {
            alignOnConfidentCmb.setSelectedIndex(0);
        } else {
            alignOnConfidentCmb.setSelectedIndex(1);
        }
    }

    /**
     * Indicates whether the user canceled the editing.
     *
     * @return a boolean indicating whether the user canceled the editing
     */
    public boolean isCanceled() {
        return canceled;
    }

    /**
     * Validates the user input.
     *
     * @return a boolean indicating whether the user input is valid
     */
    public boolean validateInput() {
        try {
            if (scoreCmb.getSelectedItem() != ModificationLocalizationScore.None) {
                Double temp = new Double(thresholdTxt.getText().trim());
                if (temp < 0 || temp > 100) {
                    JOptionPane.showMessageDialog(this, "Please verify the input for the score threshold.",
                            "Input Error", JOptionPane.ERROR_MESSAGE);
                    thresholdTxt.requestFocus();
                    return false;
                }
            }
        } catch (Exception e) {
            JOptionPane.showMessageDialog(this, "Please verify the input for the score threshold.",
                    "Input Error", JOptionPane.ERROR_MESSAGE);
            thresholdTxt.requestFocus();
            return false;
        }

        if (scoreCmb.getSelectedItem() == ModificationLocalizationScore.None) {
            int outcome = JOptionPane.showConfirmDialog(this,
                    "Disabling the probabilistic score will impair PTM localization and thus distinction\n"
                    + "between peptides. See help for more details. Continue with this setting?",
                    "Warning", JOptionPane.WARNING_MESSAGE);
            if (outcome == JOptionPane.CANCEL_OPTION || outcome == JOptionPane.CLOSED_OPTION) {
                return false;
            }
        }
        if (scoreCmb.getSelectedItem() != ModificationLocalizationScore.None && neutralLossesCmb.getSelectedIndex() == 0) {
            int outcome = JOptionPane.showConfirmDialog(this,
                    "In our experience probabilistic scores perform poorly when accounting for\n"
                    + "neutral losses. See help for more details. Continue with this setting?",
                    "Warning", JOptionPane.WARNING_MESSAGE);
            if (outcome == JOptionPane.CANCEL_OPTION || outcome == JOptionPane.CLOSED_OPTION) {
                return false;
            }
        }
        return true;
    }

    /**
     * Returns the PTM scoring preferences as set by the user.
     *
     * @return the PTM scoring preferences as set by the user
     */
    public ModificationLocalizationParameters getPtmScoringPreferences() {

        ModificationLocalizationParameters ptmScoringPreferences = new ModificationLocalizationParameters();

        ptmScoringPreferences.setProbabilisticScoreCalculation(scoreCmb.getSelectedItem() != ModificationLocalizationScore.None);
        ptmScoringPreferences.setSelectedProbabilisticScore((ModificationLocalizationScore) scoreCmb.getSelectedItem());
        ptmScoringPreferences.setProbabilisticScoreNeutralLosses(neutralLossesCmb.getSelectedIndex() == 0);
        
        if (!thresholdTxt.getText().isEmpty()) {
            ptmScoringPreferences.setProbabilisticScoreThreshold(new Double(thresholdTxt.getText().trim()));
        }

        ptmScoringPreferences.setAlignNonConfidentModifications(alignOnConfidentCmb.getSelectedIndex() == 0);

        return ptmScoringPreferences;
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jSeparator1 = new javax.swing.JSeparator();
        backgroundPanel = new javax.swing.JPanel();
        cancelButton = new javax.swing.JButton();
        okButton = new javax.swing.JButton();
        ptmScoringPanel = new javax.swing.JPanel();
        scoreCmb = new javax.swing.JComboBox();
        thresholdTxt = new javax.swing.JTextField();
        thresholdLabel = new javax.swing.JLabel();
        scoreTypeLabel = new javax.swing.JLabel();
        neutralLossesLabel = new javax.swing.JLabel();
        neutralLossesCmb = new javax.swing.JComboBox();
        siteAlignmentPanel = new javax.swing.JPanel();
        alignOnConfidentLbl = new javax.swing.JLabel();
        alignOnConfidentCmb = new javax.swing.JComboBox();
        helpJButton = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("PTM Localization");
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

        ptmScoringPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Modificatoin Scoring"));
        ptmScoringPanel.setOpaque(false);

        scoreCmb.setModel(new DefaultComboBoxModel(com.compomics.util.experiment.identification.modification.ModificationLocalizationScore.values()));
        scoreCmb.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                scoreCmbActionPerformed(evt);
            }
        });

        thresholdTxt.setHorizontalAlignment(javax.swing.JTextField.CENTER);

        thresholdLabel.setText("Threshold");

        scoreTypeLabel.setText("Probabilistic Score");

        neutralLossesLabel.setText("Account Neutral Losses");

        neutralLossesCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));

        javax.swing.GroupLayout ptmScoringPanelLayout = new javax.swing.GroupLayout(ptmScoringPanel);
        ptmScoringPanel.setLayout(ptmScoringPanelLayout);
        ptmScoringPanelLayout.setHorizontalGroup(
            ptmScoringPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(ptmScoringPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(ptmScoringPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(ptmScoringPanelLayout.createSequentialGroup()
                        .addComponent(neutralLossesLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 150, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 37, Short.MAX_VALUE)
                        .addComponent(neutralLossesCmb, javax.swing.GroupLayout.PREFERRED_SIZE, 150, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(ptmScoringPanelLayout.createSequentialGroup()
                        .addGap(0, 0, Short.MAX_VALUE)
                        .addComponent(scoreCmb, javax.swing.GroupLayout.PREFERRED_SIZE, 150, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, ptmScoringPanelLayout.createSequentialGroup()
                        .addGroup(ptmScoringPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(scoreTypeLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 150, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(thresholdLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 150, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(thresholdTxt, javax.swing.GroupLayout.PREFERRED_SIZE, 150, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap())
        );
        ptmScoringPanelLayout.setVerticalGroup(
            ptmScoringPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(ptmScoringPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(ptmScoringPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(scoreTypeLabel)
                    .addComponent(scoreCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(ptmScoringPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(neutralLossesLabel)
                    .addComponent(neutralLossesCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(ptmScoringPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(thresholdLabel)
                    .addComponent(thresholdTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap())
        );

        siteAlignmentPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Site Alignment"));
        siteAlignmentPanel.setOpaque(false);

        alignOnConfidentLbl.setText("Confident Sites");
        alignOnConfidentLbl.setToolTipText("Align peptide ambiguously localized PTMs on confident sites");

        alignOnConfidentCmb.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));

        javax.swing.GroupLayout siteAlignmentPanelLayout = new javax.swing.GroupLayout(siteAlignmentPanel);
        siteAlignmentPanel.setLayout(siteAlignmentPanelLayout);
        siteAlignmentPanelLayout.setHorizontalGroup(
            siteAlignmentPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(siteAlignmentPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(alignOnConfidentLbl)
                .addGap(18, 18, Short.MAX_VALUE)
                .addComponent(alignOnConfidentCmb, javax.swing.GroupLayout.PREFERRED_SIZE, 150, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );
        siteAlignmentPanelLayout.setVerticalGroup(
            siteAlignmentPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(siteAlignmentPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(siteAlignmentPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(alignOnConfidentLbl)
                    .addComponent(alignOnConfidentCmb, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        helpJButton.setIcon(new javax.swing.ImageIcon(getClass().getResource("/icons/help.GIF"))); // NOI18N
        helpJButton.setToolTipText("Help");
        helpJButton.setBorder(null);
        helpJButton.setBorderPainted(false);
        helpJButton.setContentAreaFilled(false);
        helpJButton.setFocusable(false);
        helpJButton.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        helpJButton.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        helpJButton.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseEntered(java.awt.event.MouseEvent evt) {
                helpJButtonMouseEntered(evt);
            }
            public void mouseExited(java.awt.event.MouseEvent evt) {
                helpJButtonMouseExited(evt);
            }
        });
        helpJButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                helpJButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout backgroundPanelLayout = new javax.swing.GroupLayout(backgroundPanel);
        backgroundPanel.setLayout(backgroundPanelLayout);
        backgroundPanelLayout.setHorizontalGroup(
            backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, backgroundPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(siteAlignmentPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, backgroundPanelLayout.createSequentialGroup()
                        .addGap(10, 10, 10)
                        .addComponent(helpJButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(okButton, javax.swing.GroupLayout.PREFERRED_SIZE, 65, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(cancelButton))
                    .addComponent(ptmScoringPanel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap())
        );
        backgroundPanelLayout.setVerticalGroup(
            backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, backgroundPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(ptmScoringPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(siteAlignmentPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addGroup(backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.CENTER)
                    .addComponent(helpJButton)
                    .addComponent(okButton)
                    .addComponent(cancelButton))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
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
     * Cancel the dialog.
     *
     * @param evt
     */
    private void cancelButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cancelButtonActionPerformed
        canceled = true;
        dispose();
    }//GEN-LAST:event_cancelButtonActionPerformed

    /**
     * Close the dialog.
     *
     * @param evt
     */
    private void okButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_okButtonActionPerformed
        if (validateInput()) {
            dispose();
        }
    }//GEN-LAST:event_okButtonActionPerformed

    /**
     * Enable/disable the subsettings.
     *
     * @param evt
     */
    private void scoreCmbActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_scoreCmbActionPerformed
        if (scoreCmb.getSelectedItem() != ModificationLocalizationScore.None) {
            neutralLossesCmb.setEnabled(true);
            thresholdTxt.setEnabled(true);
            thresholdTxt.setEditable(true);
        } else {
            neutralLossesCmb.setEnabled(false);
            thresholdTxt.setEnabled(false);
            thresholdTxt.setEditable(false);
        }
    }//GEN-LAST:event_scoreCmbActionPerformed

    /**
     * Cancel the dialog.
     *
     * @param evt
     */
    private void formWindowClosing(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowClosing
        canceled = true;
    }//GEN-LAST:event_formWindowClosing

    /**
     * Change the cursor to a hand cursor.
     *
     * @param evt
     */
    private void helpJButtonMouseEntered(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_helpJButtonMouseEntered
        setCursor(new java.awt.Cursor(java.awt.Cursor.HAND_CURSOR));
    }//GEN-LAST:event_helpJButtonMouseEntered

    /**
     * Change the cursor back to the default cursor.
     *
     * @param evt
     */
    private void helpJButtonMouseExited(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_helpJButtonMouseExited
        setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }//GEN-LAST:event_helpJButtonMouseExited

    /**
     * Open the help dialog.
     *
     * @param evt
     */
    private void helpJButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_helpJButtonActionPerformed
        setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        new HelpDialog(parentFrame, getClass().getResource("/helpFiles/PtmLocalizationPreferences.html"),
                Toolkit.getDefaultToolkit().getImage(getClass().getResource("/icons/help.GIF")),
                Toolkit.getDefaultToolkit().getImage(getClass().getResource("/icons/peptide-shaker.gif")),
                "PTM Localization - Help");
        setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }//GEN-LAST:event_helpJButtonActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JComboBox alignOnConfidentCmb;
    private javax.swing.JLabel alignOnConfidentLbl;
    private javax.swing.JPanel backgroundPanel;
    private javax.swing.JButton cancelButton;
    private javax.swing.JButton helpJButton;
    private javax.swing.JSeparator jSeparator1;
    private javax.swing.JComboBox neutralLossesCmb;
    private javax.swing.JLabel neutralLossesLabel;
    private javax.swing.JButton okButton;
    private javax.swing.JPanel ptmScoringPanel;
    private javax.swing.JComboBox scoreCmb;
    private javax.swing.JLabel scoreTypeLabel;
    private javax.swing.JPanel siteAlignmentPanel;
    private javax.swing.JLabel thresholdLabel;
    private javax.swing.JTextField thresholdTxt;
    // End of variables declaration//GEN-END:variables

}
