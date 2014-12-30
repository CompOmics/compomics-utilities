package com.compomics.software.dialogs;

import com.compomics.util.preferences.UtilitiesUserPreferences;

/**
 * Interface for parents of JavaMemoryDialog and JavaHomeDialog.
 *
 * @author Harald Barsnes
 */
public interface JavaHomeOrMemoryDialogParent {

    /**
     * Restart the given tool with the new Java options.
     */
    public void restart();

    /**
     * Returns the utilities user preferences.
     *
     * @return the utilities user preferences
     */
    public UtilitiesUserPreferences getUtilitiesUserPreferences();
}
