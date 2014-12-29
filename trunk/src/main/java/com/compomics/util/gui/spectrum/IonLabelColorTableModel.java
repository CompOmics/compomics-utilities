package com.compomics.util.gui.spectrum;

import com.compomics.util.experiment.biology.Ion;
import com.compomics.util.experiment.biology.Ion.IonType;
import com.compomics.util.experiment.biology.IonFactory;
import com.compomics.util.experiment.biology.NeutralLoss;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import javax.swing.table.DefaultTableModel;

/**
 * A table model to use for the ion label annotation colors.
 *
 * @author Marc Vaudel
 */
public class IonLabelColorTableModel extends DefaultTableModel {

    /**
     * The list of ions.
     */
    private HashMap<String, Ion> ionMap;
    /**
     * The keys.
     */
    private ArrayList<String> keys;

    /**
     * Constructor which sets a new table.
     *
     * @param iontypes the ion types
     * @param neutralLosses the neutral losses
     */
    public IonLabelColorTableModel(HashMap<IonType, HashSet<Integer>> iontypes, ArrayList<NeutralLoss> neutralLosses) {

        ionMap = new HashMap<String, Ion>();
        keys = new ArrayList<String>();
        Ion currentIon;

        for (IonType ionType : iontypes.keySet()) {
            if (ionType == IonType.IMMONIUM_ION) {
                //@TODO: we usually group immonium ions, is it a good idea?
                ionMap.put("Immonium Ion", Ion.getGenericIon(ionType, 0));
                keys.add("Immonium Ion");
            } else {
                for (Integer subtype : iontypes.get(ionType)) {
                    if (ionType == Ion.IonType.REPORTER_ION) {
                        currentIon = Ion.getGenericIon(ionType, subtype);
                        String key = currentIon.getName();
                        ionMap.put(key, currentIon);
                        keys.add(key);
                    } else {
                        for (ArrayList<NeutralLoss> possibleCombination : IonFactory.getAccountedNeutralLosses(neutralLosses)) {
                            currentIon = Ion.getGenericIon(ionType, subtype, possibleCombination);
                            String key = currentIon.getName();
                            ionMap.put(key, currentIon);
                            keys.add(key);
                        }
                    }
                }
            }
        }
    }

    /**
     * Returns the ion type at the given row.
     *
     * @param rowIndex the table row index
     * @return the ion type at the given row
     */
    public Ion getIonAtRow(int rowIndex) {
        return ionMap.get(keys.get(rowIndex));
    }

    /**
     * Constructor which sets a new empty table.
     *
     */
    public IonLabelColorTableModel() {
        ionMap = new HashMap<String, Ion>();
        keys = new ArrayList<String>();
    }

    @Override
    public int getRowCount() {
        if (keys == null) {
            return 0;
        }
        return keys.size();
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
                return "Ion";
            case 2:
                return "  ";
            default:
                return "";
        }
    }

    @Override
    public Object getValueAt(int row, int column) {
        switch (column) {
            case 0:
                return row + 1;
            case 1:
                return keys.get(row);
            case 2:
                return SpectrumPanel.determineFragmentIonColor(ionMap.get(keys.get(row)), true);
            default:
                return "";
        }
    }

    @Override
    public Class getColumnClass(int columnIndex) {
        for (int i = 0; i < getRowCount(); i++) {
            if (getValueAt(i, columnIndex) != null) {
                return getValueAt(i, columnIndex).getClass();
            }
        }
        return String.class;
    }

    @Override
    public boolean isCellEditable(int rowIndex, int columnIndex) {
        return columnIndex == 2;
    }
}
