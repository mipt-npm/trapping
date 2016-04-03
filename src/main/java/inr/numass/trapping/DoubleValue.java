/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package inr.numass.trapping;

/**
 * Класс нужен исключительно чтобы сделать простой доступ к Сишным экспортам
 * @author Darksnake
 */
public class DoubleValue {
    private double value;

    public DoubleValue(double value) {
        this.value = value;
    }

    
    
    /**
     * @return the value
     */
    public double getValue() {
        return value;
    }

    /**
     * @param value the value to set
     */
    public void setValue(double value) {
        this.value = value;
    }
}
