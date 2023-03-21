/**
 * # RNA_PECAM_DAPI

* **Developed for:** Vianney
* **Team:** Brunet
* **Date:** March 2023
* **Software:** Fiji


### Images description

3D images taken with a x63 objective on an Airyscan

2 channels:
  1. *488:* Gene1 foci
  2. *555:* vaisseaux PECAM
  3. *405:* DAPI
  4. *642:* Gene2 foci
  
A *.roi* or *.zip* file containing ROI(s) can be provided with each image.

### Plugin description

In each ROI,
* Detect Gene1 foci with Median filtering + DoG filtering + MaxEntropy thresholding
* Detect Gene2 foci with Median filtering + DoG filtering + MaxEntropy thresholding
* Detect nuclei with cellPose (cyto2)
* Find foci inside nucleus_vessel+/-

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **CellPose** 

### Version history

Version 1 released on March 21, 2023.
**/



import RNA_PECAM_DAPI_Tools.Nucleus;
import RNA_PECAM_DAPI_Tools.Tools;
import ij.*;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.ArrayList;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Objects3DIntPopulation;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;


public class RNA_PECAM_DAPI implements PlugIn {

    private final boolean canceled = false;
    public String outDirResults = "";
    private String imageDir = "";

    private RNA_PECAM_DAPI_Tools.Tools tools = new Tools();
      

    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            imageDir += IJ.getDirectory("Choose images directory")+File.separator;
            if (imageDir == null) {
                return;
            }
            File inDir = new File(imageDir);
            ArrayList<String> imageFiles = tools.findImages(imageDir, "czi");
            if (imageFiles == null) {
                return;
            }
            
            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find chanels, image calibration
            reader.setId(imageFiles.get(0));
            String[] channels = tools.findChannels(imageFiles.get(0), meta, reader);
            tools.cal = tools.findImageCalib(meta);
            String[] chs = tools.dialog(channels);
            if(chs == null)
                return;
            
            // create output folder
            outDirResults = inDir + File.separator+ "Results"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            
            // Write headers results for results files
            FileWriter fileResults = new FileWriter(outDirResults + "results.xls", false);
            BufferedWriter outPutResults = new BufferedWriter(fileResults);
            outPutResults.write("ImageName\tROI Volume\t#Nuclei\tNucleus volume (µm3)\tPECAM+\t#Foci Gene1\tFoci Gene1 volume (µm3)\t#Foci Gene2\tFoci Gene2 volume (µm3)\n");
            outPutResults.flush();            
            
            
            // Read images
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                
                // Find ROI file
                String roiFile = imageDir+rootName+".zip";
                if (!new File(roiFile).exists()) {
                    roiFile = imageDir+rootName+".roi";
                    if (!new File(roiFile).exists()) {
                        System.out.println("No ROI file found !");
                        return;
                    }
                }
               
                ImporterOptions options = new ImporterOptions();
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                options.setId(f);
                options.setSplitChannels(true);

                RoiManager rm = new RoiManager(false);
                rm.runCommand("Open", roiFile);
                Roi[] rois = rm.getRoisAsArray();

                // Open DAPI
                int indexCh = ArrayUtils.indexOf(channels, chs[2]);
                System.out.println("Opening DAPI channel = "+ chs[2]);
                ImagePlus imgDAPI = BF.openImagePlus(options)[indexCh];

                // get Object3D from rois
                Objects3DIntPopulation roisPop = tools.getObjectsFromRoi(imgDAPI, rois);
                
                // Find nucleus population
                ArrayList<Nucleus> nuclei = tools.cellposeDetection(imgDAPI, roisPop);
                tools.closeImages(imgDAPI);

                // Open gene1
                indexCh = ArrayUtils.indexOf(channels, chs[0]);
                System.out.println("Opening gene1 channel = "+ chs[0]);
                ImagePlus imgGene1 = BF.openImagePlus(options)[indexCh];
                Objects3DIntPopulation gene1Pop = tools.findGenesPop(imgGene1, nuclei, 1);
                System.out.println(gene1Pop.getNbObjects()+" genes1 found in nucleus");
                tools.closeImages(imgGene1);

                // Open gene2
                indexCh = ArrayUtils.indexOf(channels, chs[1]);
                System.out.println("Opening gene2 channel = "+ chs[1]);
                ImagePlus imgGene2 = BF.openImagePlus(options)[indexCh];
                Objects3DIntPopulation gene2Pop = tools.findGenesPop(imgGene2, nuclei, 2);
                System.out.println(gene2Pop.getNbObjects()+" genes2 found in nucleus");

                // Write parameters
                IJ.showStatus("Writing parameters ...");
                double roisVol = tools.findPopVolume(roisPop);
                for (Nucleus nucleus: nuclei) {
                    outPutResults.write(rootName+"\t"+roisVol+"\t"+nucleus.params.get("index")+"\t"+nucleus.params.get("nucVol")+"\t"+nucleus.params.get("vessel")
                        +"\t"+nucleus.params.get("fociGene1")+"\t"+nucleus.params.get("fociGene1Vol")
                        +"\t"+nucleus.params.get("fociGene2")+"\t"+nucleus.params.get("fociGene2Vol")+"\n");
                outPutResults.flush();
                }

                // save image objects
                IJ.showStatus("Save images objects ...");
                String path = outDirResults + rootName+"_Objects.tif";
                tools.saveGenesImage(nuclei, gene1Pop, gene2Pop, roisPop, imgGene2, path);  
                imgGene2.close();
        }
            outPutResults.close();
            IJ.showStatus("Process done");
        } catch (DependencyException | ServiceException | FormatException | IOException ex) {
            Logger.getLogger(RNA_PECAM_DAPI.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
