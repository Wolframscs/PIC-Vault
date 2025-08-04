package macro;
import java.util.*;
import star.common.*;
import star.base.neo.*;
import star.vis.*;
import star.meshing.*;
import star.base.report.*;
import java.io.File;

public class batchMeshInterface extends StarMacro{

    public void execute()
    {
    Simulation simulation_0=getActiveSimulation();
    
    simulation_0.get(MeshOperationManager.class).executeAll();
    simulation_0.getRegionManager().updateInterfacesFromPartContacts(new NeoObjectVector(new Object[] {simulation_0.get(SimulationPartManager.class)}), RegionManager.CreateInterfaceMode.CONTACT);
// AutoMesh   
    String originalPath = simulation_0 .getSessionPath();
    File file = new File(originalPath);
    String baseName = file.getName().replace(".sim", "");
    String newName = baseName + "_mesh.sim";
    String newPath = file.getParent() + File.separator + newName;
    simulation_0.saveState(newPath);
// Model Runing  
    simulation_0.getSimulationIterator().run();
//    String originalPath = simulation_0 .getSessionPath();
//    File file = new File(originalPath);
//    String baseName = file.getName().replace(".sim", "");
    String newName2 = baseName + "_result.sim";
    String newPath2 = file.getParent() + File.separator + newName2;
    simulation_0.saveState(newPath2);
//  simulation_0.saveState(resolvePath("./NBC17_TR_Model_NoBatSim_NoMesh_Save_result.sim"));
    }
}
