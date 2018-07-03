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
	
	var files = Directory.GetFiles(@"E:\databases\PDB\cif\complete");
	var locker = new object(); //pro paralelní běh
	var counter = 0; // počítadlo souborů, které už jsem prošla
	
    HashSet<String> ligandy = new HashSet<String>();
	ligandy = File.ReadAllLines(@"vyhledavani_ligandu/vyrazene_ligandy.txt").ToHashSet(); //načítání souboru s ligandy
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
		var str = StructureReader.Read(e).Structure;  // nacti strukturu
	
		var listOfRes = str.PdbResidues().Select(a => a.Name).Distinct().ToList(); //výběr residuí z načtené struktury

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
			if ((counter % 2) == 0){ //ať nevypisuju po jednom (zase jen pro přehled během běhu programu)
				counter.Dump();
			}
		}
	}catch(Exception ex) {
		ex.Dump();
	}
	});
	
	File.WriteAllLines(@"pdb_copy_files2.txt", copyFiles);//neprojíté struktury
	File.WriteAllLines(@"pdb_structures_zbytek.txt", structures);//vybrané struktury
}