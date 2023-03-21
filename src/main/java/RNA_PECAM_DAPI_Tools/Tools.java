package RNA_PECAM_DAPI_Tools;




import RNA_PECAM_DAPI.Cellpose.CellposeSegmentImgPlusAdvanced;
import RNA_PECAM_DAPI.Cellpose.CellposeTaskSettings;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackMerge;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Point3D;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.VoxelInt;
import mcib3d.geom2.measurements.MeasureCentroid;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import org.apache.commons.io.FilenameUtils;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;



/**
 *
 * @author phm
 */

public class Tools {
    
    public CLIJ2 clij2 = CLIJ2.getInstance();
    
    // min size for dots
    private double minFoci = 0.05;
    private double maxFoci = 50;
    private final double minDOGFoci = 1;
    private final double maxDOGFoci = 2;
    private String geneThreshold = "MaxEntropy";
    
    
     // Cellpose
    public String cellPoseModel = "cyto2";
    public String cellPoseEnvDirPath = (IJ.isWindows()) ? System.getProperty("user.home")+"\\miniconda3\\envs\\CellPose" : "/opt/miniconda3/envs/cellpose";
    double minNucVol = 20;
    double maxNucVol = Double.MAX_VALUE;
    
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    public Calibration cal = new Calibration();    
    public float pixVol = 0;
    
     /**
     * check  installed modules
     * @return 
     */
    public boolean checkInstalledModules() {
        // check install
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.log("CLIJ not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }

        
   /**
     * Find images in folder
     */
    public ArrayList findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No Image found in "+imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }
    
    
    /**
     * Dialog
     */
    public String[] dialog(String[] channels) {
        String[] channelsName = {"Gene1 : ", "Gene2 : ", "DAPI : ", "Vessel : "};
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 20, 0);
        gd.addImage(icon);
        gd.addMessage("Channels selection", Font.getFont("Monospace"), Color.blue);
        for (int n = 0; n < channelsName.length; n++) {
            gd.addChoice(channelsName[n], channels, channels[0]);
        }
        gd.addMessage("Dots filter", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min foci volume : ", minFoci, 2, 6, "µm3");
        gd.addNumericField("Max foci volume : ", maxFoci, 2, 6, "µm3");
      
        // Calibration
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY pixel size : ", cal.pixelWidth, 3);
        gd.addNumericField("Z pixel size : ", cal.pixelDepth, 3);
        gd.showDialog();
        String[] ch = new String[channelsName.length];
        for (int i = 0; i < channelsName.length; i++)
            ch[i] = gd.getNextChoice();

        minFoci = gd.getNextNumber();
        maxFoci= gd.getNextNumber();
        cal.pixelWidth = gd.getNextNumber();
        cal.pixelDepth = gd.getNextNumber();        
        pixVol = (float) (cal.pixelWidth*cal.pixelHeight*cal.pixelDepth);
        if (gd.wasCanceled())
                ch = null;
        return(ch);
    } 
    
    /**
     * Find channels name
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs+1];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n).toString();
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelFluor(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelFluor(0, n).toString();
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break;    
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break; 
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        channels[chs] = "None";
        return(channels);     
    }
    
    /**
     * Find image calibration
     * @param meta
     * @return 
     */
    public Calibration findImageCalib(IMetadata meta) {
        // read image calibration
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);
        return(cal);
    }
    
    
    /**
     *
     * @param img
     */
    public void closeImages(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    
     /**
     * Difference of Gaussians 
     * Using CLIJ2
     * @param img
     * @param size1
     * @param size2
     * @return imgGauss
     */ 
    public ImagePlus DOG(ImagePlus img, double size1, double size2) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLDOG = clij2.create(imgCL);
        clij2.differenceOfGaussian3D(imgCL, imgCLDOG, size1, size1, size1, size2, size2, size2);
        clij2.release(imgCL);
        ImagePlus imgDOG = clij2.pull(imgCLDOG);
        clij2.release(imgCLDOG);
        return(imgDOG);
    }
    
     /**
     * Median filter using CLIJ2
     * @param img
     * @param sizeXY
     * @param sizeZ
     * @return 
     */ 
    public ImagePlus median_filter(ImagePlus img, double sizeXY) {
       ClearCLBuffer imgCL = clij2.push(img); 
       ClearCLBuffer imgCLMed = clij2.create(imgCL);
       clij2.median3DBox(imgCL, imgCLMed, sizeXY, sizeXY, sizeXY);
       clij2.release(imgCL);
       ImagePlus imgMed = clij2.pull(imgCLMed);
       clij2.release(imgCLMed);
       return(imgMed);
    } 
    
    /**
     * Threshold 
     * USING CLIJ2
     * @param img
     * @param thMed
     * @return 
     */
    public ImagePlus threshold(ImagePlus img, String thMed) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        clij2.release(imgCL);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        clij2.release(imgCLBin);
        return(imgBin);
    }
    
    /**
     * return objects population in an binary image
     * @param img
     * @return pop
     */
    public Objects3DIntPopulation getPopFromImage(ImagePlus img, Calibration cal) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        labels.setCalibration(cal);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(labels);
        return pop;
    }
    
    
    
     /**
     * Remove object with size < min and size > max
     * @param pop
     * @param min
     * @param max
     */
    public void popFilterSize(Objects3DIntPopulation pop, double min, double max) {
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < min) || (new MeasureVolume(p).getVolumeUnit() > max));
        pop.resetLabels();
    }
    
     /**
     * Clear out side roi
     * @param img
     * @param roi
     */
    public void clearOutSide(ImagePlus img, Roi roi) {
        PolygonRoi poly = new PolygonRoi(roi.getFloatPolygon(), Roi.FREEROI);
        poly.setLocation(0, 0);
        for (int n = 1; n <= img.getNSlices(); n++) {
            ImageProcessor ip = img.getImageStack().getProcessor(n);
            ip.setRoi(poly);
            ip.setBackgroundValue(0);
            ip.fillOutside(poly);
        }
        img.deleteRoi();
        img.updateAndDraw();
    }
    
    /**
     * Return Object3D from roi
     */
    public Objects3DIntPopulation getObjectsFromRoi(ImagePlus img, Roi[] rois) {
        ImagePlus imgDup = new Duplicator().run(img);
        Objects3DIntPopulation roisPop = new Objects3DIntPopulation();
        for (Roi roi : rois) {
            imgDup.setRoi(roi);
            IJ.run(imgDup, "8-bit","");
            IJ.setBackgroundColor(0, 0, 0);
            IJ.run(imgDup, "Clear Outside", "stack");
            IJ.setForegroundColor(255, 255, 255);
            IJ.run(imgDup, "Fill", "stack");
            imgDup.deleteRoi();
            ImageHandler imh = ImageHandler.wrap(imgDup);
            Object3DInt obj = new Object3DInt(imh, 255);
            roisPop.addObject(obj);
        }
        closeImages(imgDup);
        return(roisPop);  
    }
    
    
     /**
     * Look for all 2D cells: 
     * - apply CellPose in 2D slice 
     * 
     * @param img
     * @return 
     * @throws java.io.IOException
     */
    public ArrayList<Nucleus> cellposeDetection(ImagePlus img, Objects3DIntPopulation roisPop) throws IOException{
        ArrayList<Nucleus> nuclei = new ArrayList<>();
        ImagePlus imgResized = img.resize((int)(img.getWidth()*0.5), (int)(img.getHeight()*0.5), 1, "none");
        ImagePlus imgMed = median_filter(imgResized, 1);
        closeImages(imgResized);
        // Define CellPose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(cellPoseModel, 1, 10, cellPoseEnvDirPath);
        settings.useGpu(true);
        settings.setStitchThreshold(0.5);
       
        // Run CellPose
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgMed);
        ImagePlus imgOut = cellpose.run();
        imgOut = imgOut.resize(img.getWidth(), img.getHeight(), "none");
        imgOut.setCalibration(cal);
        closeImages(imgMed);
        // Get cells as a population of objects
        Objects3DIntPopulation nucPop = new Objects3DIntPopulation(ImageHandler.wrap(imgOut));
        System.out.println(nucPop.getNbObjects() + " cell detections");
        popFilterSize(nucPop, minNucVol, maxNucVol);
        System.out.println(nucPop.getNbObjects() + " cells remaining after size filtering"); 
        
        // tag nucleus inside vessel (roi)
        for (Object3DInt nucObj : nucPop.getObjects3DInt()) {
            VoxelInt center = new MeasureCentroid(nucObj).getCentroidRoundedAsVoxelInt();
            Nucleus nucleus = new Nucleus(nucObj);
            nuclei.add(nucleus);
            nucleus.params.put("index", (double)nucObj.getLabel());
            nucleus.params.put("nucVol", new MeasureVolume(nucObj).getVolumeUnit());
            for (Object3DInt roiObj : roisPop.getObjects3DInt()) {
                if (roiObj.contains(center)) {
                    nucleus.params.put("vessel", 1.0);
                }
                else
                    nucleus.params.put("vessel", 0.0);
                break;
            }
        }
        closeImages(imgOut);
        return(nuclei);
    } 
    
    /**
     * Find coloc between pop1 and pop2
     * tag nucleus with foci number and foci volume 
     * @param nuclei
     * @param fociPop
     * @throws java.io.IOException
     */
    public Objects3DIntPopulation findColocPop (ArrayList<Nucleus> nuclei, Objects3DIntPopulation fociPop, int gene) throws IOException {
        Objects3DIntPopulation fociNucPop = new Objects3DIntPopulation();
        float fociIndex = 0;
        for (Nucleus nucleus: nuclei) {
            Object3DInt nucObj = nucleus.nucleus;
            double foci = 0;
            double fociVol = 0;
            for (Object3DInt fociObj : fociPop.getObjects3DInt()) {
                Point3D center = new MeasureCentroid(fociObj).getCentroidAsPoint();
                if (nucObj.contains(new VoxelInt(center.getRoundX(), center.getRoundY(), center.getRoundZ(), 255))) {
                    foci++;
                    fociVol += new MeasureVolume(fociObj).getVolumeUnit();
                    Object3DInt fociNuc = fociObj;
                    fociIndex++;
                    fociNuc.setLabel(fociIndex);
                    fociNucPop.addObject(fociNuc);
                }
            }
            switch (gene) {
                case 1 :
                    nucleus.params.put("fociGene1", foci);
                    nucleus.params.put("fociGene1Vol", fociVol);
                    break;
                case 2 :
                    nucleus.params.put("fociGene2", foci);
                    nucleus.params.put("fociGene2Vol", fociVol);
                    break;
            }
        }
        fociNucPop.resetLabels();
        return(fociNucPop);
    }

    
    /**
     * Find genes population
     * @param imgGene
     * @return genePop
     */
    public Objects3DIntPopulation findGenesPop(ImagePlus imgGene, ArrayList<Nucleus> nuclei, int gene) throws IOException {
        IJ.showStatus("Finding gene dots ...");
        ImagePlus imgDup = new Duplicator().run(imgGene);
        ImagePlus imgDOG = DOG(imgDup, minDOGFoci, maxDOGFoci);
        closeImages(imgDup);
        ImagePlus imgBin = threshold(imgDOG, geneThreshold); 
        closeImages(imgDOG);
        imgBin.setCalibration(cal);
        Objects3DIntPopulation genePop = getPopFromImage(imgBin, cal);
        popFilterSize(genePop, minFoci, maxFoci);
        System.out.println(genePop.getNbObjects() + " genes"+gene+" found");
        closeImages(imgBin);
        // tag nuclei with dots number and volume
        Objects3DIntPopulation colocPop = findColocPop(nuclei, genePop, gene);
        return(colocPop);
    }
    

    /**
     * Find sum volume of objects  
     * @param dotsPop
     * @return vol
     */
    
    public double findPopVolume (Objects3DIntPopulation dotsPop) {
        IJ.showStatus("Findind object's volume");
        double sumVol = 0;
        for(Object3DInt obj : dotsPop.getObjects3DInt()) {
            sumVol += new MeasureVolume(obj).getVolumeUnit();
        }
        return(sumVol);
    }
    
    
   /**
     * Find roi volume
     * @param roi
     * @param img
     * @return volume
     */
    public double roiVolume(Roi roi, ImagePlus imgAstro) {
        PolygonRoi poly = new PolygonRoi(roi.getFloatPolygon(), PolygonRoi.FREEROI); 
        poly.setLocation(0, 0);
        ImageProcessor ip = imgAstro.getProcessor();
        ip.setRoi(poly);
        ImageStatistics stats = ImageStatistics.getStatistics(ip, Measurements.AREA, cal);
        double volume = stats.area * imgAstro.getNSlices();
        return(volume);
    }
    
    
    /**
     * save images objects population
     * @param img
     * @param gene1Pop
     * @param gene2Pop
     * @param path
     */
    public void saveGenesImage (ArrayList<Nucleus> nuclei, Objects3DIntPopulation gene1Pop, Objects3DIntPopulation gene2Pop, Objects3DIntPopulation roisPop, 
            ImagePlus img, String path) {
        // red gene1 , green gene2, blue nuclei, grey roi
        ImageHandler imhDapi = ImageHandler.wrap(img).createSameDimensions();
        ImageHandler imhGene1 = imhDapi.createSameDimensions();
        ImageHandler imhGene2 = imhDapi.createSameDimensions();
        ImageHandler imhRoi = imhDapi.createSameDimensions();
        
        // draw nuclei
        for (Nucleus nuc : nuclei)
            nuc.nucleus.drawObject(imhDapi, 255);
        // draw genes population
        for (Object3DInt ob : gene1Pop.getObjects3DInt())
            ob.drawObject(imhGene1, 255);
        for (Object3DInt ob : gene2Pop.getObjects3DInt())
            ob.drawObject(imhGene2, 255);
        for (Object3DInt ob : roisPop.getObjects3DInt())
            ob.drawObject(imhRoi, 32); 
        ImagePlus[] imgColors = {imhGene1.getImagePlus(), imhGene2.getImagePlus(), imhDapi.getImagePlus(), imhRoi.getImagePlus()};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, true);
        imgObjects.setCalibration(cal);
        // Save images
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(path);
        imhGene1.closeImagePlus();
        imhGene2.closeImagePlus();
        imhDapi.closeImagePlus();
        imhRoi.closeImagePlus();
    }
   
}
