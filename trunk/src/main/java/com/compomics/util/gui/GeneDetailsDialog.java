package com.compomics.util.gui;

import com.compomics.util.examples.BareBonesBrowserLaunch;
import com.compomics.util.experiment.annotation.gene.GeneFactory;
import com.compomics.util.experiment.annotation.go.GOFactory;
import com.compomics.util.experiment.identification.matches.ProteinMatch;
import java.awt.Component;
import java.awt.event.MouseEvent;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import javax.swing.JPanel;
import javax.swing.JTable;
import javax.swing.ScrollPaneConstants;
import javax.swing.border.TitledBorder;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import no.uib.jsparklines.extra.HtmlLinksRenderer;

/**
 * This dialog displays the gene details associated to a protein match.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class GeneDetailsDialog extends javax.swing.JDialog {

    /**
     * The GO terms factory.
     */
    private GOFactory goFactory = GOFactory.getInstance();
    /**
     * The Gene factory.
     */
    private GeneFactory geneFactory = GeneFactory.getInstance();
    /**
     * the protein accessions of this match.
     */
    private ArrayList<String> proteinAccessions;
    /**
     * The GO term descriptions attached to this protein match.
     */
    private ArrayList<String> goTermDescriptions;
    /**
     * The color to use for the HTML tags for the selected rows, in HTML color
     * code.
     */
    private String selectedRowHtmlTagFontColor = "#FFFFFF"; // @TODO: move somewhere more generic...
    /**
     * The color to use for the HTML tags for the rows that are not selected, in
     * HTML color code.
     */
    private String notSelectedRowHtmlTagFontColor = "#0101DF"; // @TODO: move somewhere more generic...

    /**
     * Creates a new GeneDetailsDialog.
     *
     * @param parent the parent frame
     * @param proteinMatchKey the protein match key
     * @throws IOException
     */
    public GeneDetailsDialog(java.awt.Frame parent, String proteinMatchKey) throws IOException {
        super(parent, true);
        initComponents();
        proteinAccessions = new ArrayList<String>(Arrays.asList(ProteinMatch.getAccessions(proteinMatchKey)));
        goTermDescriptions = goFactory.getProteinGoDescriptions(proteinMatchKey);
        Collections.sort(goTermDescriptions);
        setUpGUI();
        setLocationRelativeTo(parent);
        setVisible(true);
    }

    /**
     * Set up the GUI.
     */
    private void setUpGUI() {

        goTable.getColumn("Accession").setCellRenderer(new HtmlLinksRenderer(selectedRowHtmlTagFontColor, notSelectedRowHtmlTagFontColor));

        // set the preferred size of the accession column
        Integer width = getPreferredAccessionColumnWidth(goTable, goTable.getColumn("Accession").getModelIndex(), 20);
        if (width != null) {
            goTable.getColumn("Accession").setMinWidth(width);
            goTable.getColumn("Accession").setMaxWidth(width);
        } else {
            goTable.getColumn("Accession").setMinWidth(15);
            goTable.getColumn("Accession").setMaxWidth(Integer.MAX_VALUE);
        }

        goTable.getColumn(" ").setMaxWidth(50);
        goTable.getColumn(" ").setMinWidth(50);

        goTable.getTableHeader().setReorderingAllowed(false);
        goTable.setAutoCreateRowSorter(true);

        // correct the color for the upper right corner
        JPanel proteinCorner = new JPanel();
        proteinCorner.setBackground(goTable.getTableHeader().getBackground());
        goTableScrollPane.setCorner(ScrollPaneConstants.UPPER_RIGHT_CORNER, proteinCorner);

        // make sure that the scroll panes are see-through
        goTableScrollPane.getViewport().setOpaque(false);

        try {
            String title = "", geneIdsTxt = "", geneNamesTxt = "", chromosomeTxt = "";
            ArrayList<String> geneNames = new ArrayList<String>();
            for (String accession : proteinAccessions) {
                if (title.equals("")) {
                    title += "Gene details for ";
                } else {
                    title += ", ";
                }
                title += accession;

                String geneName = geneFactory.getGeneNameForUniProtProtein(accession);
                geneNames.add(geneName);
            }

            ArrayList<String> chromosomes = new ArrayList<String>();
            for (String geneName : geneNames) {
                if (!geneIdsTxt.equals("")) {
                    geneIdsTxt += ", ";
                    geneNamesTxt += ", ";
                }

                String ensemblId = geneFactory.getGeneEnsemblId(geneName);

                if (ensemblId == null) {
                    geneIdsTxt += "unknown";
                } else {
                    geneIdsTxt += geneFactory.getGeneEnsemblId(geneName);
                }

                if (geneName == null) {
                    geneNamesTxt += "unknown";
                } else {
                    geneNamesTxt += geneName;
                }

                String chromosome = geneFactory.getChromosomeForGeneName(geneName);
                chromosomes.add(chromosome);
            }

            for (String chromosome : chromosomes) {
                if (!chromosomeTxt.equals("")) {
                    chromosomeTxt += ", ";
                }
                if (chromosome == null) {
                    chromosomeTxt += "unknown";
                } else {
                    chromosomeTxt += chromosome;
                }
            }

            ((TitledBorder) detailsPanel.getBorder()).setTitle(title);
            geneIdTxt.setText(geneIdsTxt);
            geneNameTxt.setText(geneNamesTxt);
            chromosomeNameTxt.setText(chromosomeTxt);

        } catch (ClassNotFoundException e) { // @TODO: better error handling
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (IllegalArgumentException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
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
        detailsPanel = new javax.swing.JPanel();
        ensemlbIdLabel = new javax.swing.JLabel();
        geneIdTxt = new javax.swing.JTextField();
        geneNameTxt = new javax.swing.JTextField();
        geneNameLabel = new javax.swing.JLabel();
        chromosomeLabel = new javax.swing.JLabel();
        chromosomeNameTxt = new javax.swing.JTextField();
        goAnnotationLabel = new javax.swing.JLabel();
        goTableScrollPane = new javax.swing.JScrollPane();
        goTable = new javax.swing.JTable();
        okButton = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Gene Details");

        backgroundPanel.setBackground(new java.awt.Color(230, 230, 230));

        detailsPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Gene Details"));
        detailsPanel.setOpaque(false);

        ensemlbIdLabel.setText("Ensembl Gene ID");

        geneIdTxt.setEditable(false);

        geneNameTxt.setEditable(false);

        geneNameLabel.setText("Gene Name");

        chromosomeLabel.setText("Chromosome");

        chromosomeNameTxt.setEditable(false);

        goAnnotationLabel.setText("GO Annotation");

        goTable.setModel(new GOTableModel());
        goTable.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseReleased(java.awt.event.MouseEvent evt) {
                goTableMouseReleased(evt);
            }
            public void mouseExited(java.awt.event.MouseEvent evt) {
                goTableMouseExited(evt);
            }
        });
        goTable.addMouseMotionListener(new java.awt.event.MouseMotionAdapter() {
            public void mouseMoved(java.awt.event.MouseEvent evt) {
                goTableMouseMoved(evt);
            }
        });
        goTableScrollPane.setViewportView(goTable);

        javax.swing.GroupLayout detailsPanelLayout = new javax.swing.GroupLayout(detailsPanel);
        detailsPanel.setLayout(detailsPanelLayout);
        detailsPanelLayout.setHorizontalGroup(
            detailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(detailsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(detailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(goTableScrollPane, javax.swing.GroupLayout.DEFAULT_SIZE, 707, Short.MAX_VALUE)
                    .addGroup(detailsPanelLayout.createSequentialGroup()
                        .addGroup(detailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(ensemlbIdLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 110, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(geneNameLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 110, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(detailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(geneNameTxt)
                            .addComponent(geneIdTxt)))
                    .addGroup(detailsPanelLayout.createSequentialGroup()
                        .addComponent(goAnnotationLabel)
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addGroup(detailsPanelLayout.createSequentialGroup()
                        .addComponent(chromosomeLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 110, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(chromosomeNameTxt)))
                .addContainerGap())
        );
        detailsPanelLayout.setVerticalGroup(
            detailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(detailsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(detailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(ensemlbIdLabel)
                    .addComponent(geneIdTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(detailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(geneNameLabel)
                    .addComponent(geneNameTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(detailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(chromosomeLabel)
                    .addComponent(chromosomeNameTxt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addComponent(goAnnotationLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(goTableScrollPane, javax.swing.GroupLayout.DEFAULT_SIZE, 206, Short.MAX_VALUE)
                .addContainerGap())
        );

        okButton.setText("OK");
        okButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                okButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout backgroundPanelLayout = new javax.swing.GroupLayout(backgroundPanel);
        backgroundPanel.setLayout(backgroundPanelLayout);
        backgroundPanelLayout.setHorizontalGroup(
            backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(backgroundPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(detailsPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, backgroundPanelLayout.createSequentialGroup()
                        .addGap(0, 0, Short.MAX_VALUE)
                        .addComponent(okButton, javax.swing.GroupLayout.PREFERRED_SIZE, 86, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap())
        );
        backgroundPanelLayout.setVerticalGroup(
            backgroundPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, backgroundPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(detailsPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(okButton)
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
     * Close the dialog.
     *
     * @param evt
     */
    private void okButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_okButtonActionPerformed
        dispose();
    }//GEN-LAST:event_okButtonActionPerformed

    /**
     * If the user clicks the GO accession column the GO term is opened in the
     * web browser.
     *
     * @param evt
     */
    private void goTableMouseReleased(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_goTableMouseReleased
        int row = goTable.getSelectedRow();
        int column = goTable.getSelectedColumn();

        if (row != -1) {
            this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));

            if (evt == null || evt.getButton() == MouseEvent.BUTTON1) {

                // open protein link in web browser
                if (column == goTable.getColumn("Accession").getModelIndex() && evt.getButton() == MouseEvent.BUTTON1
                        && ((String) goTable.getValueAt(row, column)).lastIndexOf("<html>") != -1) {

                    String link = (String) goTable.getValueAt(row, column);
                    link = link.substring(link.indexOf("\"") + 1);
                    link = link.substring(0, link.indexOf("\""));

                    this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
                    BareBonesBrowserLaunch.openURL(link);
                    this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
                }
            }

            this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
        }
    }//GEN-LAST:event_goTableMouseReleased

    /**
     * Changes the cursor into a hand cursor if the table cell contains an HTML
     * link.
     *
     * @param evt
     */
    private void goTableMouseMoved(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_goTableMouseMoved
        int row = goTable.rowAtPoint(evt.getPoint());
        int column = goTable.columnAtPoint(evt.getPoint());

        if (column == goTable.getColumn("Accession").getModelIndex() && goTable.getValueAt(row, column) != null) {

            String tempValue = (String) goTable.getValueAt(row, column);

            if (tempValue.lastIndexOf("<html>") != -1) {
                this.setCursor(new java.awt.Cursor(java.awt.Cursor.HAND_CURSOR));
            } else {
                this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
            }
        } else {
            this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
        }
    }//GEN-LAST:event_goTableMouseMoved

    /**
     * Changes the cursor back to the default cursor a hand.
     *
     * @param evt
     */
    private void goTableMouseExited(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_goTableMouseExited
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }//GEN-LAST:event_goTableMouseExited
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JPanel backgroundPanel;
    private javax.swing.JLabel chromosomeLabel;
    private javax.swing.JTextField chromosomeNameTxt;
    private javax.swing.JPanel detailsPanel;
    private javax.swing.JLabel ensemlbIdLabel;
    private javax.swing.JTextField geneIdTxt;
    private javax.swing.JLabel geneNameLabel;
    private javax.swing.JTextField geneNameTxt;
    private javax.swing.JLabel goAnnotationLabel;
    private javax.swing.JTable goTable;
    private javax.swing.JScrollPane goTableScrollPane;
    private javax.swing.JButton okButton;
    // End of variables declaration//GEN-END:variables

    /**
     * Table model for the GO annotation.
     */
    private class GOTableModel extends DefaultTableModel {

        @Override
        public int getRowCount() {
            if (goTermDescriptions != null) {
                return goTermDescriptions.size();
            }
            return 0;
        }

        @Override
        public int getColumnCount() {
            return 3;
        }

        @Override
        public String getColumnName(int column) {

            switch (column) {
                case 0:
                    return " ";
                case 1:
                    return "Accession";
                case 2:
                    return "Description";
                default:
                    return "";
            }
        }

        @Override
        public Object getValueAt(int row, int column) {

            switch (column) {
                case 0:
                    return (row + 1);
                case 1:
                    try {
                        return addGoLink(goFactory.getTermAccession(goTermDescriptions.get(row)));
                    } catch (Exception e) {
                        return "Error";
                    }
                case 2:
                    try {
                        String accession = goFactory.getTermAccession(goTermDescriptions.get(row));
                        return goFactory.getTermDescription(accession);
                    } catch (Exception e) {
                        return "Error";
                    }
                default:
                    return "";
            }
        }

        @Override
        public Class getColumnClass(int columnIndex) {
            return getValueAt(0, columnIndex).getClass();
        }

        @Override
        public boolean isCellEditable(int rowIndex, int columnIndex) {
            return false;
        }
    }

    /**
     * Returns the GO accession number as a web link to the given GO term at
     * QuickGO.
     *
     * @param goAccession
     * @return the GO accession number as a web link to the given GO term at
     * QuickGO
     */
    public String addGoLink(String goAccession) { // @TODO: move method to utilities...
        return "<html><a href=\"" + getGoAccessionLink(goAccession)
                + "\"><font color=\"" + notSelectedRowHtmlTagFontColor + "\">"
                + goAccession + "</font></a></html>";
    }

    /**
     * Returns the GO accession number as a web link to the given GO term at
     * QuickGO.
     *
     * @param goAccession the GO accession number
     * @return the GO accession web link
     */
    public String getGoAccessionLink(String goAccession) {
        return "http://www.ebi.ac.uk/QuickGO/GTerm?id=" + goAccession;
    }

    /**
     * Gets the preferred width of the column specified by colIndex. The column
     * will be just wide enough to show the column head and the widest cell in
     * the column. Margin pixels are added to the left and right (resulting in
     * an additional width of 2*margin pixels. Returns null if the max width
     * cannot be set.
     *
     * @param table the table
     * @param colIndex the colum index
     * @param margin the margin to add
     * @return the preferred width of the column
     */
    public Integer getPreferredAccessionColumnWidth(JTable table, int colIndex, int margin) {

        DefaultTableColumnModel colModel = (DefaultTableColumnModel) table.getColumnModel();
        TableColumn col = colModel.getColumn(colIndex);

        // get width of column header
        TableCellRenderer renderer = col.getHeaderRenderer();
        if (renderer == null) {
            renderer = table.getTableHeader().getDefaultRenderer();
        }

        Component comp = renderer.getTableCellRendererComponent(table, col.getHeaderValue(), false, false, 0, 0);
        int width = comp.getPreferredSize().width;

        // add margin
        width += 2 * margin;

        return width;
    }
}
