<Query Kind="Expression">
  <Reference>E:\Zuzka\CH_stacking\Newtonsoft.Json.dll</Reference>
  <Reference>E:\Zuzka\CH_stacking\System.Interactive.dll</Reference>
  <Reference>&lt;RuntimeDirectory&gt;\System.IO.dll</Reference>
  <Reference>E:\Zuzka\CH_stacking\System.Reactive.Core.dll</Reference>
  <Reference>&lt;RuntimeDirectory&gt;\System.Threading.dll</Reference>
  <Reference>&lt;RuntimeDirectory&gt;\System.Threading.Tasks.dll</Reference>
  <Reference>&lt;RuntimeDirectory&gt;\System.Threading.Tasks.Parallel.dll</Reference>
  <Reference>&lt;RuntimeDirectory&gt;\System.Threading.Timer.dll</Reference>
  <Reference>E:\Zuzka\CH_stacking\WebChemistry.Framework.Core.dll</Reference>
  <Reference>E:\Zuzka\CH_stacking\WebChemistry.Queries.Core.dll</Reference>
  <Namespace>System.Threading.Tasks</Namespace>
  <Namespace>WebChemistry.Framework.Core</Namespace>
  <Namespace>WebChemistry.Framework.Core.Pdb</Namespace>
  <Namespace>WebChemistry.Framework.Math</Namespace>
  <Namespace>WebChemistry.Queries.Core</Namespace>
  <IncludePredicateBuilder>true</IncludePredicateBuilder>
</Query>

using System.Threading.Tasks;
using WebChemistry.Framework.Core;
using WebChemistry.Queries.Core;
using WebChemistry.Queries.Service;
using System.Activities.Expressions;

/*
    rozlišení glykosilace ... druhé kolo filtrování
 */
void Main()
{
	Directory.SetCurrentDirectory(@"E:\Zuzka\CH_stacking");//výběr pracovního adresáře
	
	//použít jen jedno z nich!!!!
	var files = Directory.GetFiles("motivy/C3C4/HIS/all/motives/"); //seznam PDB motivů nalezených při úvodním prohledávání
	var structures = File.ReadAllLines(@"motivy/C3C4/HIS/pdb_HIS_all.csv").ToHashSet(); //vstupní soubor obsahuje cesty k pdb strukturám ... použít jen vyfiltrovaný seznam pdb
	structures.Count.Dump();//test
	
	//zápis do souborů
	var zapisLigandy = new HashSet<String>();
	var zapisGl1 = new HashSet<String>();
	var zapisGlykosilace = new HashSet<String>();
	var pdbGlykosilace = new HashSet<String>();
	var pdbSpornaGlykosilace = new HashSet<String>();
	var pdbLigandy = new HashSet<String>();
	var pdbGl1 = new HashSet<String>();
	
    /*
        definice patternQuery dotazů
     */
	var tetrose = QueryBuilder.Rings(new string[]{"C", "C", "C", "C", "O"});
	var pentose = QueryBuilder.Rings(new string[]{"C", "C", "C", "C", "C", "O"});
													
	//dotaz pro získání motivu použitý v původním vyhledávání (všechny struktury by měly mít nějaký motiv)
	var query01 = QueryBuilder.
			Or(tetrose, pentose).
			Spherify(5).ConnectedResidues(0).
			Filter(l => QueryBuilder.Count(l, QueryBuilder.Residues(new string[]{"HIS"})) > 0);
			
	var query = query01.ToMetaQuery().Compile();
	
	var sugar = QueryBuilder.Or(tetrose, pentose).ToMetaQuery().Compile();
	
	//jednoduchá glykosilace (dotaz, jestli je obsahuje cukr, který je navázaný na protein)
	var connection = QueryBuilder.Or(tetrose, pentose).AmbientResidues(1.5).
					 Filter (l => QueryBuilder.Count(l, QueryBuilder.AminoAcids()) > 0).ToMetaQuery().Compile();			
	
    //sleduje navázané residua k cukru v blízkosti aromatické aminokyseliny
	var test = QueryBuilder.Path(query01, QueryBuilder.HetResidues(), QueryBuilder.AminoAcids());
	
    /*
        hlavičky pro .csv soubory
     */
	zapisLigandy.Add("pdb_id;name;motive_signatures;");
	zapisGl1.Add("pdb_id;name;motive_signatures;");
	zapisGlykosilace.Add("pdb_id;motive_name;lenght;motive_signatures;");
	pdbGlykosilace.Add("pdb_id;count;");
	pdbSpornaGlykosilace.Add("pdb_id;count;");

    /*
        počítadla
     */
	var countLigands = 0;
	var countGl1 = 0;
	var totalCounter = 0;
	var gl = false; //testovací proměnná pro rozlišení glykosilace
	
    //procházení struktur
	foreach (var e in structures)
	{
	try{
		var str = StructureReader.Read(e).Structure;  // nacti strukturu
		var result = query.Matches(str); //porovnání s query
		var motivesCount = 0;//počítadlo nalezených motivů
		var name = new StringBuilder(Path.GetFileNameWithoutExtension(e));//jméno pdb pro zápis motivu

        /*
            formátování zápisu do souboru
         */
		var writeLigands = new StringBuilder();
		var writeGl1 = new StringBuilder();
		var writeGlycosilation = new StringBuilder();
		
        //procházení nalezených motivů ve struktuře
		foreach (var r in result){
			name.Append(motivesCount);//přidání čísla motivu do názvu !!!v aktuální verzi to druhý motiv zapíše jako 01, třetí 012 atd.!!!
			var str2 = r.ToStructure(name.ToString(), true, true); // převedení motivu do pdb struktury
			
            //test, jestli struktura motivu obsahuje cukr přímo navázaný na protein
            var gl1 = connection.Matches(str2);
			
			if (gl1.Count > 0){//na cukr je přímo navázaná AMK
				/*
                    zápis .pdb motivu do kategorie jednoduché glykosilace
                    zápis textu + struktury motivu
                */
				writeGl1.Append(str.Id);
				writeGl1.Append(";");
				writeGl1.Append(name.ToString());
				writeGl1.Append(";");
				writeGl1.Append(r.Signature);
				writeGl1.Append(";");
				var pathGl1 = new StringBuilder("motivy/C3C4/HIS/all/gl1/");
				pathGl1.Append(str.Id);
				pathGl1.Append(motivesCount);
				using (TextWriter write = File.CreateText(pathGl1.ToString())){
					str2.WritePdb(write);
				}
				pdbGl1.Add(str.Id);
				zapisGl1.Add(writeGl1.ToString());
				countGl1++;
			}else{ 
				var sugars = sugar.Matches(str2); //najde cukry ve struktuře motivu
                /*
                    hledám residua navázaná na cukru, postupně rozšiřuju patternQuery dotaz
                    !!! neuvažuju případ, že by v daném motivu bylo víc cukrů - testuju jen ten první!!!
                 */
				var deepSearch = new StringBuilder(); //formulace patternQuery dotazu
				var resID = sugars.First().Atoms.First().ResidueIdentifier();
				deepSearch.Append("Path(ResidueIds(\"");
				deepSearch.Append(resID);
				deepSearch.Append("\"), HetResidues()");
				
                /*
                    postupné rozšiřování dotazu, dokud je v řetězci něco navázáno a není to aminokyselina 
                    100 je rezerva, které by to nikdy nemělo dosáhnout
                */
				for (int i = 1; i < 100; i++){
					var countOfGlycosilation = 0;
					var testBuilderAA = new StringBuilder();//test, jestli v aktuálním řetězci existuje spojení s AMK
					testBuilderAA.Append(deepSearch.ToString());
					testBuilderAA.Append(", AminoAcids())");
					
					var testQueryAA = PythonEngine.GetQuery(testBuilderAA.ToString());
					var glTestAA = testQueryAA.Matches(str);//testuju, jestli je motiv v PŮVODNÍ struktuře (motiv může být delší než okolí 5A) 
                                                                //=> teoreticky to může vytvořit duplicitu
					
					var testBuilder = new StringBuilder();//test, jestli existuje oligosacharid
					testBuilder.Append(deepSearch.ToString());
					testBuilder.Append(")");
					var testQuery = PythonEngine.GetQuery(testBuilder.ToString());
					var glTest = testQuery.Matches(str); //znovu testuju v PŮVODNÍ struktuře (ne  ve struktuře nalezeného motivu)
					
					if (glTest.Count == 0 && i == 1){//na cukr není navázaný žádný další cukr (dál už nehledám)
						/*
                            zápis .pdb motivu do kategorie ligandy
                            zápis textu + struktury motivu
                        */
						writeLigands.Append(str.Id);
						writeLigands.Append(";");
						writeLigands.Append(name.ToString());
						writeLigands.Append(";");
						writeLigands.Append(r.Signature);
						writeLigands.Append(";");
						var pathLigands = new StringBuilder("motivy/C3C4/HIS/all/ligand/");
						pathLigands.Append(str.Id);
						pathLigands.Append(motivesCount);
						using (TextWriter write = File.CreateText(pathLigands.ToString())){
							str2.WritePdb(write);
						}
						countLigands++;
						pdbLigandy.Add(str.Id);
						zapisLigandy.Add(writeLigands.ToString());
						gl = false; //rozlišení, jestli je motiv v kateforii glykosylace
						break;
					}else if (glTest.Count > 0 && glTestAA.Count == 0){//existuje oligosacharid, ale ještě jsme nenašli, jestli je navázaný na AMK NEBO je to jen ligand!!!
						deepSearch.Append(", HetResidues()");
						
						var connection2 = new StringBuilder(); //test, jestli tam opravdu není vazba na AMK (test pomocí AmbientResidues)
						
                        //procházení motivů, které obsahují řetězec HetResiduí začínající cukrem
						foreach (var motive in glTest){
                            /*
                                test přítomnosti vazby pomocí Ambient Residues (důsledek častých chyb v PDB databázi)
                             */
							var rID = motive.Atoms.Last().ResidueIdentifier();//ID posledního residua v motivu
							connection2.Append("ResidueIds(\"");
							connection2.Append(rID);
							connection2.Append("\").AmbientResidues(1.5).Filter(lambda m: m.Count(AminoAcids())>0)");
							var connectionQuery = PythonEngine.GetQuery(connection2.ToString());
							var cQ = connectionQuery.Matches(str);
							if (cQ.Count > 0){ //motiv měl být v kategorii glykosylace, ale nenašla jsem ho kvůli chybě v PDB databáze
                                /*
                                    zápis .pdb motivu do kategorie glykosylace
                                    JEN zápis textu, protože pokračuju v hledání
                                 */
								writeGlycosilation.Append(str.Id);
								writeGlycosilation.Append(";");
								writeGlycosilation.Append(name.ToString());
								writeGlycosilation.Append(";");
								writeGlycosilation.Append(i+1);
								writeGlycosilation.Append(";");
								writeGlycosilation.Append(r.Signature);
								writeGlycosilation.Append(";");
								countOfGlycosilation++;
								zapisGlykosilace.Add(writeGlycosilation.ToString());
								gl = true;
							}
							connection2 = new StringBuilder();
						}
					}else if (glTest.Count > 0 && glTestAA.Count > 0){//našla jsem glykosylaci oligosacharidu, ale hledám dál, jestli nenajdu větší
                        /*
                            zápis .pdb motivu do kategorie glykosylace
                            JEN zápis textu, protože pokračuju v hledání
                        */
						writeGlycosilation.Append(str.Id);
						writeGlycosilation.Append(";");
						writeGlycosilation.Append(name.ToString());
						writeGlycosilation.Append(";");
						writeGlycosilation.Append(i+1);
						writeGlycosilation.Append(";");
						writeGlycosilation.Append(r.Signature);
						writeGlycosilation.Append(";");
						deepSearch.Append(", HetResidues()");
						countOfGlycosilation++;
						zapisGlykosilace.Add(writeGlycosilation.ToString());
						gl = true;
						
					}else{//už to nic nenašlo => předchozí výsledek byl maximální
						if(gl){//rozlišení ligandu od glykosylace
                            /*
                                zápis .pdb motivu do kategorie glykosylace
                                zápis textu + struktury
                            */
							var pathGlycosilation2 = new StringBuilder("motivy/C3C4/HIS/all/glykosylace/");
							pathGlycosilation2.Append(str.Id);
							pathGlycosilation2.Append(motivesCount);
							using (TextWriter write = File.CreateText(pathGlycosilation2.ToString())){
								str2.WritePdb(write);
							}
							var zapis = new StringBuilder();
							zapis.Append(str.Id);
							zapis.Append(";");
							zapis.Append(result.Count);
							zapis.Append(";");
							pdbGlykosilace.Add(zapis.ToString());
						}else{
							/*
                                zápis .pdb motivu do kategorie ligandy
                                zápis textu + struktury
                            */
							writeLigands.Append(str.Id);
							writeLigands.Append(";");
							writeLigands.Append(name.ToString());
							writeLigands.Append(";");
							writeLigands.Append(r.Signature);
							writeLigands.Append(";");
							var pathLigands2 = new StringBuilder("motivy/C3C4/HIS/all/ligand/");
							pathLigands2.Append(str.Id);
							pathLigands2.Append(motivesCount);
							using (TextWriter write = File.CreateText(pathLigands2.ToString())){
								str2.WritePdb(write);
							}
							countLigands++;
							pdbLigandy.Add(str.Id);
							zapisLigandy.Add(writeLigands.ToString());
						}
						gl = false;
						break;
					    }
				    }
			    }
			    motivesCount++;
		    }
		    totalCounter++;
		    totalCounter.Dump();
	    }catch(Exception ex){
		    ex.Dump();
	    }
	}
	
	"počet lidandů".Dump();
	countLigands.Dump();
	"počet gl1".Dump();
	countGl1.Dump();
	
	File.WriteAllLines(@"HIS_ligands.csv", zapisLigandy);
	File.WriteAllLines(@"HIS_gl1.csv", zapisGl1);
	File.WriteAllLines(@"HIS_glycosilation.csv", zapisGlykosilace);
}
