/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.compomics.util.experiment.identification.protein_inference.fm_index;

import com.compomics.util.experiment.personalization.ExperimentObject;

/**
 *
 * @author dominik.kopczynski
 */
public class SNPElement extends ExperimentObject  {
    public int position;
    public char sourceAA;
    public char targetAA;
    

    
    public SNPElement(int _position, char source, char target){
        position = _position;
        sourceAA = source;
        targetAA = target;
    }
}
