//doplnit
using WebChemistry.Queries.Core;
using WebChemistry.Framework.Core;
using System.Threading.Tasks;

/*
    výběr pdb struktur s daným ligandem
 */

void Main()
{
	Directory.SetCurrentDirectory(@"E:\Zuzka\CH_stacking"); // vyber pracovni adresar
	
	var structures = new List<string>();//budu tam ukládat vybrané pdb struktury
	
	//var files = Directory.GetFiles(@"Z:\999986-WebChemData\Databases\PDB\Assemblies");
	var files = Directory.GetFiles(@"Z:\999986-WebChemData\Databases\PDB\mirrors\units", "*.*", SearchOption.AllDirectories);

	//procházení doplňkových stuktur ze souboru
	/*HashSet<String> problemove_struktury_puvodni = File.ReadAllLines(@"E:\Zuzka\CH_stacking\vyber_pdb\16_05_2019\pdb_copy_files.txt").ToHashSet();
	HashSet<String> problemove_struktury = new HashSet<String>();
	foreach (var i in problemove_struktury_puvodni){
		var kod = i.Split('\\');
		//kod[7].Substring(0,4).Dump();
		problemove_struktury.Add(kod[7].Substring(0,4));
	}*/
	
	var locker = new object(); //pro paralelní běh
	var counter = 0; // počítadlo souborů, které už jsem prošla
	
    HashSet<String> ligandy = new HashSet<String>();
	ligandy = File.ReadAllLines(@"vyhledavani_ligandu/C3C4/ligandsC3C4.txt").ToHashSet(); //načítání souboru s ligandy
	var ligandName = new HashSet<String>();//soubor s ligandy neobsahuje jen název ligandu => úprava
	
    /*
        obejítí OutOfMemoryException
        možná už nebude potřeba
    */
	var copyFiles = new HashSet<String>();
	foreach (var e in files){
		copyFiles.Add(e);
	}
	
    /*
        úprava textu ze souboru s ligandy, abych dostala jen název ligandu
     */
	foreach (var x in ligandy){
		var pole = x.Split('_');
		ligandName.Add(pole[0]);
	}
	
	//projítí pdb databáze
	Parallel.ForEach(files, new ParallelOptions(){MaxDegreeOfParallelism = 4}, e =>
	{try{
		//projdu jen struktury ze souboru copy_files
		/*var pdb_kod = e.Split('\\');
		pdb_kod[5].Substring(0,4).Dump();
		HashSet<String> ans = new HashSet<String>(problemove_struktury);
		HashSet<String> pdb = new HashSet<String>();
		pdb.Add(pdb_kod[5].Substring(0,4));
		ans.IntersectWith(pdb);
		if (ans.Count() != 0){*/

		var str = StructureReader.Read(e).Structure;  // nacti strukturu
	
		var listOfRes = str.PdbResidues().Select(a => a.Name).Distinct().ToList(); //výběr residuí z načtené struktury

	    str = null;
        /*
            porovním seznam vybraných residuí se seznamem residuí aktuálně načtené pdb struktury
        */	
		var prunik = listOfRes.Intersect(ligandName); 
		if (prunik.Count() != 0){
			e.Dump();//pro přehled v běhu programu
			structures.Add(e);
		}
		
		lock (locker){
			copyFiles.Remove(e);//pojistka, kdyby to spadlo, tak vím, které struktury jsem prošla
			counter++;
			if ((counter % 100000) == 0){ //ať nevypisuju po jednom (zase jen pro přehled během běhu programu)
				counter.Dump();

				//Z paměťových důvodů mi to padá dřív, než se dostanu k zápisu struktur, takže zkusím zapisovat průběžně
				var structuresFile = new StringBuilder("pdb_structures_");
				structuresFile.Append(counter);
				structuresFile.Append(".txt");
				File.WriteAllLines(structuresFile.ToString(), structures);
				
				var copyFile = new StringBuilder("pdb_copy_files_");
				copyFile.Append(counter);
				copyFile.Append(".txt");
				File.WriteAllLines(copyFile.ToString(), copyFiles);
			}
		}
	}catch(Exception ex) {
		ex.Dump();
	}
	});
	
	File.WriteAllLines(@"pdb_copy_files.txt", copyFiles);//neprojíté struktury
	File.WriteAllLines(@"pdb_structures.txt", structures);//vybrané struktury
}