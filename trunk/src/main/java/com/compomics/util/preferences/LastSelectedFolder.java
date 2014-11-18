/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.compomics.util.preferences;

import java.io.Serializable;
import java.util.HashMap;

/**
 * Convenience class keeping class of the last selected folders
 *
 * @author Marc
 */
public class LastSelectedFolder implements Serializable {

    /**
     * Serial version UID for backward compatibility.
     */
    static final long serialVersionUID = 7612497861284156966L;
    /**
     * Map of the last selected folders indexed by use cases.
     */
    private HashMap<String, String> lastSelectedFolder = null;
    /**
     * The default use case
     */
    private static final String defaultUseCase = "default";

    /**
     * Constructs a new last selected folder class. Defaults to the user home by
     * default.
     */
    public LastSelectedFolder() {
        setLastSelectedFolder("user.home");
    }

    /**
     * Sets the last selected folder of the given use case. Default use case is
     * used if null.
     *
     * @param useCase the use case, can be null
     * @param folderPath the path to the last selected folder
     */
    public void setLastSelectedFolder(String useCase, String folderPath) {
        if (folderPath != null) {
            if (useCase == null) {
                useCase = defaultUseCase;
            }
            if (lastSelectedFolder == null) {
                lastSelectedFolder = new HashMap<String, String>(1);
            }
            lastSelectedFolder.put(useCase, folderPath);
        }
    }

    /**
     * Sets the last selected folder.
     *
     * @param folderPath the path to the last selected folder
     */
    public void setLastSelectedFolder(String folderPath) {
        setLastSelectedFolder(defaultUseCase, folderPath);
    }

    /**
     * Returns the last selected folder according to the given use case. The
     * default use case is used if usecase is null.
     *
     * @param useCase the use case
     *
     * @return the last selected folder, null if not found
     */
    public String getLastSelectedFolder(String useCase) {
        if (lastSelectedFolder == null) {
            return null;
        }
        if (useCase == null) {
            useCase = defaultUseCase;
        }
        return lastSelectedFolder.get(useCase);
    }

    /**
     * Returns the last selected folder according to the given use case.
     *
     * @return the last selected folder, null if not found
     */
    public String getLastSelectedFolder() {
        return getLastSelectedFolder(null);
    }

}
