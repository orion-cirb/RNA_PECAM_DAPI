package RNA_PECAM_DAPI_Tools;

import java.util.HashMap;
import mcib3d.geom2.Object3DInt;

/**
 * @author hm
 */
public class Nucleus {
    
    public Object3DInt nucleus;
    public HashMap<String, Double> params;
    
    public Nucleus(Object3DInt nucleus) {
        this.nucleus = nucleus;
        this.params = new HashMap<>();
    }
    
    public void setParams(double index, double nucVol, double vessel, double fociGene1, double fociGene1Vol, double fociGene2, double fociGene2Vol) {
        params.put("index", index);
        params.put("nucVol", nucVol);
        params.put("vessel", vessel);
        params.put("fociGene1", fociGene1);
        params.put("fociGene2Vol", fociGene2Vol);
        params.put("fociGene2", fociGene2);
        params.put("fociGene2Vol", fociGene2Vol);
    }
}
