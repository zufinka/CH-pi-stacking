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
	
    HashSet<String> structures = new HashSet<String>();
	structures = File.ReadAllLines(@"vyber_pdb/pdb_structures_C3C4.txt").ToHashSet(); //vstupní soubor obsahuje cesty k pdb strukturám
	structures.Count.Dump();//test
	
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

	foreach (var e in structures){
	    try{
		    var str = StructureReader.Read(e).Structure;  // nacti strukturu
		    var result = query.Matches(str); //porovnání s query
		    var zapis = new StringBuilder();//formátování výstupu do souboru
		
		    if (result.Count != 0){
			    pdb_list.Add(e);
			    zapis.Append(str.Id + ";");
			    zapis.Append(result.Count);
			    zapis.Append(";");
			    var name = new StringBuilder(Path.GetFileNameWithoutExtension(e));//jméno pdb pro zápis motivu
			    var counter = 0; //počítadlo motivů v jedné struktuře
				
				foreach (var m in result){//projdu motivy ve struktuře
					var ligandPath = new StringBuilder("motivy/C3C4/HIS/all/motives/");//formátování názvu souboru, do kterého uložím nalezené motivy
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
	    }catch(Exception ex){
		    ex.Dump();
	    }
	}
	
	File.WriteAllLines(@"pdb_ligands_HIS_all.csv", zapisLigandy);
	File.WriteAllLines(@"pdb_HIS_all.csv", pdb_list);
}
