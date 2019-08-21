using System;
using System.Threading.Tasks;
using WebChemistry.Framework.Core;
using WebChemistry.Queries.Core;

/*
    Projde vybrané pdb struktury a zjistí, jestli obsahují cukr v blízkosti aromatické aminokyseliny.
    Potřeba upravit pro každou aminokyselinu.
 */
void Main()
{
Directory.SetCurrentDirectory(@"E:\Zuzka\CH_stacking");//nastaví pracovní adresář
	var pdb_list = new HashSet<String>();//seznam odkazů na struktury, které vyhovují kritériím
	
    /*HashSet<String> structures = new HashSet<String>();
	structures = File.ReadAllLines(@"E:\Zuzka\CH_stacking\vyber_pdb\16_05_2019\final_pdb_structures.txt").ToHashSet(); //vstupní soubor obsahuje cesty k pdb strukturám
	structures.Count.Dump();//test*/

	HashSet<String> structures_upravene = new HashSet<String>();

	//vytvoří novou adresu pro přístup k souborům
	//999986-WebChemData/Databases/PDB/bio/ - úložiště pro biological assembly
	/*foreach (var i in structures){
		var deleni = i.Split('\\');
		if (deleni.Length > 7){
			String a = String.Format("{0}\\{1}\\{2}\\{3}\\bio\\{4}",deleni[0],deleni[1],deleni[2],deleni[3],deleni[7].Substring(0,8));
			structures_upravene.Add(a);
		}else{
			String a = String.Format("{0}\\{1}\\{2}\\{3}\\bio\\{4}.cif",deleni[0],deleni[1],deleni[2],deleni[3],deleni[5].Substring(0,4));
			structures_upravene.Add(a);
		}
	}*/
	
	structures_upravene = File.ReadAllLines(@"E:\Zuzka\CH_stacking\motivy\16_05_2019\HIS\copystructures4.txt").ToHashSet();
	structures_upravene.Count.Dump();
	
	var copyFiles = new HashSet<String>();
	foreach (var e in structures_upravene){
		copyFiles.Add(e);
	}
	
    /*
        patternQuery dotaz, jestli struktura obsahuje cukrový kruh a uvedenou aromatickou aminokyselinu do vzdálenosti 5A
     */
	var tetrose = QueryBuilder.Rings(new string[]{"C", "C", "C", "C", "O"});
	var pentose = QueryBuilder.Rings(new string[]{"C", "C", "C", "C", "C", "O"});
	
	var query = QueryBuilder.
			Or(tetrose, pentose).
			Spherify(5).ConnectedResidues(0).
			Filter(l => QueryBuilder.Count(l, QueryBuilder.Residues(new string[]{"HIS"})) > 0).ToMetaQuery().Compile(); 
	
    /*
        proměnná pro formátování výstupu do souboru
     */
    var zapisLigandy = new HashSet<String>();
	zapisLigandy.Add("pdb_id;motive_count;ligands_count;motive_signatures;");
	IStructure str = null;
	var pojistka = 0;

	foreach (var e in structures_upravene){
	    try{
		    str = StructureReader.Read(e).Structure;  // nacti strukturu
		    var result = query.Matches(str); //porovnání s query
		    var zapis = new StringBuilder();//formátování výstupu do souboru
			str.Id.Dump();
		    if (result.Count != 0){
			    pdb_list.Add(e);
			    zapis.Append(str.Id + ";");
			    zapis.Append(result.Count);
			    zapis.Append(";");
			    var name = new StringBuilder(Path.GetFileNameWithoutExtension(e));//jméno pdb pro zápis motivu
			    var counter = 0; //počítadlo motivů v jedné struktuře
				
				foreach (var m in result){//projdu motivy ve struktuře
					var ligandPath = new StringBuilder(@"E:\Zuzka\CH_stacking\motivy\16_05_2019\HIS\motivy\");//formátování názvu souboru, do kterého uložím nalezené motivy
					ligandPath.Append(str.Id);
					ligandPath.Append(counter);
					zapis.Append(m.Signature);
					zapis.Append(";");
				
                    /*
                        zápis nalezeného motivu do .pdb souboru
                     */
				  	using (TextWriter write = File.CreateText(ligandPath.ToString())){
						name.Append(counter);
						m.ToStructure(name.ToString(), true, true).WritePdb(write);
					}
					counter++;
				}
				zapisLigandy.Add(zapis.ToString());
		    }	
			str = null;
			copyFiles.Remove(e);
			pojistka++;
			if (pojistka % 4 == 0){
				String b = String.Format("copystructures{0}.txt", pojistka);
				File.WriteAllLines(b, copyFiles);
			}
	    }catch(Exception ex){
		    ex.Dump();
			str = null;
	    }
	}
	
	File.WriteAllLines(@"pdb_ligands_HIS_all.csv", zapisLigandy);
	File.WriteAllLines(@"pdb_HIS_all.csv", pdb_list);
}