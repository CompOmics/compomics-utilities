package com.compomics.util.general;

import com.compomics.util.experiment.personalization.ExperimentObject;

/**
 * Convenience class allowing the boxing of an object.
 *
 * @author Marc Vaudel
 *
 * @param <K> the type of object to box
 */
public class BoxedObject<K> extends ExperimentObject {

    /**
     * The object to box.
     */
    private K object;

    /**
     * Constructor.
     */
    public BoxedObject() {

    }
    /**
     * Constructor.
     * 
     * @param initialValue the initial value
     */
    public BoxedObject(K initialValue) {
        object = initialValue;
    }

    /**
     * Returns the boxed object.
     *
     * @return the boxed object
     */
    public K getObject() {
        
        return object;
    }

    /**
     * Sets the boxed object.
     *
     * @param object the boxed object
     */
    public void setObject(K object) {
        
        this.object = object;
    }
}
