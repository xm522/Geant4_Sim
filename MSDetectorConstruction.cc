
#include "MSDetectorConstruction.hh"


//#include "OpNoviceDetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"


//---------------------------------------------------------------------
/*
OpNoviceDetectorConstruction::OpNoviceDetecorConstrution()
 : G4VUserDetectorConstruction()
 {
	 ExpHall_x = ExpHall_y = ExpHall_z = 10.0*m;
	 tank_x = tank_y = tank_z = 5.0*m;
 }
 */
 MSDetectorConstruction::MSDetectorConstruction()
{
  expHall_x = expHall_y = expHall_z = 10.0*m;
  tank_x    = tank_y    = tank_z    =  5.0*m;
  bubble_x  = bubble_y  = bubble_z  =  0.5*m;
}
 //---------------------------------------------------------------------
 /*OpNoviceDetecorConstrution::~OpNoviceDetecorConstrution(){;}*/
 MSDetectorConstruction::~MSDetectorConstruction(){;}
 
//----------------------------------------------------------------------
/*G4VPhysicalVolume* OpNoviceDetecorConstrution::Construct()*/
G4VPhysicalVolume* MSDetectorConstruction::Construct()
{
	//----- Materials ----------
	
	G4double a, z, density;
	G4int nelements;
	
	// Air
	
	G4Element* N = new G4Element("Nitrogen", "N" , z = 7 , a = 14.01*g/mole);
	G4Element* O = new G4Element("Oxygen"  , "O" , z = 8 , a = 16.00*g/mole);
	
	G4Material* Air = new G4Material("Air", density = 1.29*mg/cm3, nelements=2);
	Air->AddElement(N, 70.*perCent);
	Air->AddElement(O, 30.*perCent);

	
	// Polystyrene
	
	G4Element* H = new G4Element("Hydrogen", "H", z=1, a = 1.01*g/mole);
	G4Element* C = new G4Element("Carbon"  , "C", z=6, a = 12.01*g/mole);
	
	G4Material* polystyrene = new G4Material("Polystyrene", density = 1.05*g/cm3, nelements=2);
	polystyrene->AddElement(H, 8);
	polystyrene->AddElement(C, 8);
	
	// Acrylic cladding
	G4Material* Acrylic = new G4Material("Acrylic", density = 1.20*g/cm3, nelements=3);
	Acrylic->AddElement(H, 8);
	Acrylic->AddElement(C, 5); 
	Acrylic->AddElement(O, 2);
	
//------- Generate & Add Material Properties Table -------------------------------------

const G4int nEntries = 32;

G4double photonEnergy[nEntries] =
            { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
              2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
              2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
              2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
              2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
              3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
              3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
              3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV 
             };

//---Acrylic

G4double RefractiveIndex3[nEntries]=
			{1.49, 1.49,  1.49, 1.49,  1.49,
              1.49,  1.49, 1.49,  1.49, 1.49,
              1.49, 1.49, 1.49,   1.49, 1.49,
              1.49, 1.49, 1.49, 1.49, 1.49,
              1.49, 1.49,  1.49, 1.49,  1.49,
              1.49, 1.49,  1.49, 1.49,  1.49,
              1.49,   1.49
             };
G4MaterialPropertiesTable* myMPT3 = new G4MaterialPropertiesTable(); 
  myMPT3->AddProperty("RINDEX", photonEnergy, RefractiveIndex3, nEntries);
  
  Acrylic->SetMaterialPropertiesTable(myMPT3);

//---Polystyrene

G4double RefractiveIndex1[nEntries]=
			{1.5435, 1.544,  1.5445, 1.545,  1.5455,
              1.546,  1.5465, 1.547,  1.5475, 1.548,
              1.5485, 1.5492, 1.55,   1.5505, 1.551,
              1.5518, 1.5522, 1.5530, 1.5535, 1.554,
              1.5545, 1.555,  1.5555, 1.5560,  1.5565,
              1.557, 1.5575,  1.5580, 1.5585,  1.5990,
              1.5995,   1.60
             };

G4double absorption[nEntries] =
          {13.448*mm,  14.082*mm,  16.329*mm,  19.174*mm, 12.346*mm, 13.889*mm,
           15.152*mm, 17.241*mm, 18.868*mm, 20.000*mm, 26.316*mm, 35.714*mm,
           45.455*mm, 47.619*mm, 52.632*mm, 52.632*mm, 55.556*mm, 52.632*mm,
           52.632*mm, 47.619*mm, 45.455*mm, 41.667*mm, 37.037*mm, 33.333*mm,
           30.000*mm, 28.500*mm, 27.000*mm, 24.500*mm, 22.000*mm, 19.500*mm,
           17.500*mm, 14.500*mm
          };
          

G4double scintilFast[nEntries] =
            { 2.70, 2.70, 2.70, 2.70, 2.70, 2.70, 2.70,
              2.70, 2.70, 2.70, 2.70, 2.70, 2.70, 2.70,
              2.70, 2.70, 2.70, 2.70, 2.70, 2.70, 2.70,
              2.70, 2.70, 2.70, 2.70, 2.70, 2.70, 2.70,
              2.70, 2.70, 2.70, 2.70 };
              
	
//Air - more
  G4double RefractiveIndex2[nEntries] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };

  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", photonEnergy, RefractiveIndex2, nEntries);
  
  Air->SetMaterialPropertiesTable(myMPT2);
              

G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

myMPT1->AddProperty("RINDEX",        photonEnergy, RefractiveIndex1, nEntries)
      ->SetSpline(true);
myMPT1->AddProperty("ABSLENGTH",     photonEnergy, absorption,       nEntries)
	  ->SetSpline(true);
myMPT1->AddProperty("FASTCOMPONENT", photonEnergy,scintilFast,       nEntries)
	  ->SetSpline(true);
	  
myMPT1->AddConstProperty("SCINTILLATIONYIELD",50./MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT1->AddConstProperty("FASTTIMECONSTANT", 2.0*ns);
  myMPT1->AddConstProperty("SLOWTIMECONSTANT",10.0*ns);
  myMPT1->AddConstProperty("YIELDRATIO",0.8);

polystyrene->SetMaterialPropertiesTable(myMPT1);

// ---- my volumes and stuff :D ----

// hall

/*G4Box* expHall_box = new G4Box("world",expHall_x,expHall_y,expHall_z);
 
G4LogicalVolume* expHall_logic
	= new G4LogicalVolume(expHall_box,0,"world",0,0,0);
	
G4VPhysicalVolume* expHall_phys
	= new G4PVPlacement(0,G4ThreeVector(),expHall_logic,"World",0,false,0);
*/

 G4Box* expHall_box = new G4Box("World",expHall_x,expHall_y,expHall_z);

  G4LogicalVolume* expHall_log
    = new G4LogicalVolume(expHall_box,Air,"World",0,0,0);

  G4VPhysicalVolume* expHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);



// Transformations 
G4ThreeVector posV1= G4ThreeVector(0,0,0);
 G4RotationMatrix rotm = G4RotationMatrix();
 rotm.rotateX(90.*deg);
 rotm.rotateY(0.*deg);
 rotm.rotateZ(0.*deg);
 G4Transform3D transform1 = G4Transform3D(rotm,posV1);

// cladding 

G4double pRmin2 = 0.5*mm, pRmax2 = 0.505*mm;
G4double pDz   = 5.5*cm;
	      G4double pSPhi1 = 0.0*deg, pDPhi1 = 360.0*deg;
	      G4double pSPhi2 = 0.0*deg, pDPhi2 = 360.0*deg;
G4Tubs* Cladding_V =
	      new G4Tubs("V Cladding",  //Its name
		     pRmin2, pRmax2, pDz, pSPhi1, pDPhi1);

		 
      G4LogicalVolume* CladV =
	      new G4LogicalVolume(Cladding_V, //the solid
				  Acrylic, //material
				  "V_Clad");
      
      G4VPhysicalVolume* CladV_Phys = new G4PVPlacement(transform1,                       //no rotation
                    CladV,             //its logical volume
                    "V_Clad",                //its name
                    expHall_log,                //its mother  volume
                    false,                   //yes boolean operation
                    0,                       //copy number
                    0);   

// fibers

G4double pRmin = 0.*mm, pRmax = 0.5*mm;
//	      G4double pDz   = 5.5*cm;
//	      G4double pSPhi1 = 0.0*deg, pDPhi1 = 360.0*deg;
//	      G4double pSPhi2 = 0.0*deg, pDPhi2 = 360.0*deg;
      G4Tubs* WireV =
	      new G4Tubs("V Wire",  //Its   name
		     pRmin, pRmax, pDz, pSPhi1, pDPhi1);

      G4LogicalVolume* logicWireV =
	      new G4LogicalVolume(WireV, //the solid
				  polystyrene, //material
				  "V_Wire");
      
      G4VPhysicalVolume* WireV_Phys = new G4PVPlacement(transform1,                       //no rotation
                    logicWireV,             //its logical volume
                    "V_Wire",                //its name
                    expHall_log,                //its mother  volume
                    false,                   //yes boolean operation
                    0,                       //copy number
                    0);   

//------Mesh Wire Horizontal -------------------------------------------------------
 //G4NistManager*nist = G4NistManager::Instance();
 //G4Material*WireH_mat= nist->FindOrBuildMaterial("G4_POLYSTYRENE");
 G4ThreeVector posH1 = G4ThreeVector(-1.01*mm,-1.01*mm,0);
	
	G4Tubs* WireH =
	      new G4Tubs("H Wire",  //Its name
	      pRmin, pRmax, pDz, pSPhi2, pDPhi2); 

 G4LogicalVolume* logicWireH =
	      new G4LogicalVolume(WireH, //the solid
				  polystyrene, //material
				  "H_Wire");
      
    G4VPhysicalVolume* WireH_Phys =   new G4PVPlacement(0,                       //no rotation
                    posH1,                    //at position
                    logicWireH,             //its logical volume
                    "H_Wire",                //its name
                    expHall_log,                //its mother  volume
                    false,                   //yes boolean operation
                    0,                       //copy number
                    0);   
//horizontal cladding
  
G4Tubs* Cladding_H =
	      new G4Tubs("H Cladding",  //Its name
		     pRmin2, pRmax2, pDz, pSPhi1, pDPhi1);

		 
      G4LogicalVolume* CladH =
	      new G4LogicalVolume(Cladding_V, //the solid
				  Acrylic, //material
				  "V_Clad");
      
      G4VPhysicalVolume* CladH_Phys = new G4PVPlacement(0,posH1,                       //no rotation
                    CladH,             //its logical volume
                    "V_Clad",                //its name
                    expHall_log,                //its mother  volume
                    false,                   //yes boolean operation
                    0,                       //copy number
                    0);   
// surfaces 

G4OpticalSurface* OplogicWireSurface = new G4OpticalSurface("WireSurface");
  OplogicWireSurface->SetType(dielectric_dielectric);
  OplogicWireSurface->SetFinish(polished);
  OplogicWireSurface->SetModel(unified);
  
new G4LogicalBorderSurface("WireSurfaceV", 
			WireV_Phys,expHall_phys, OplogicWireSurface);
			
new G4LogicalBorderSurface("WireSurfaceH", 
			WireH_Phys,expHall_phys, OplogicWireSurface);
			
new G4LogicalBorderSurface("CladSurfaceV", 
			CladV_Phys,expHall_phys, OplogicWireSurface);
new G4LogicalBorderSurface("CladSurfaceH", 
			CladH_Phys,expHall_phys, OplogicWireSurface);
const G4int num = 2;
  G4double Ephoton[num] = {2.034*eV, 4.136*eV};

  //OplogicWireSurface
  G4double RefractiveIndex[num] = {1.35, 1.40};
  G4double SpecularLobe[num]    = {0.3, 0.3};
  G4double SpecularSpike[num]   = {0.2, 0.2};
  G4double Backscatter[num]     = {0.2, 0.2};
			
G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();
  
  myST1->AddProperty("RINDEX",                Ephoton, RefractiveIndex, num);
  myST1->AddProperty("SPECULARLOBECONSTANT",  Ephoton, SpecularLobe,    num);
  myST1->AddProperty("SPECULARSPIKECONSTANT", Ephoton, SpecularSpike,   num);
  myST1->AddProperty("BACKSCATTERCONSTANT",   Ephoton, Backscatter,     num);

  OplogicWireSurface->SetMaterialPropertiesTable(myST1);

return expHall_phys;
}
