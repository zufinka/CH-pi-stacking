// potřeba načíst namespaces, bez toho to nebude fungovat
using System.Threading.Tasks;
using WebChemistry.Framework.Core;
using WebChemistry.Queries.Core;
/*
    výběr ligandů, které obsahují pentosovy a hexosovy cyklus bez násobné vazby
*/

void Main(){	
	/*
    	výběr pracovního adresáře
 	*/
	Directory.SetCurrentDirectory(@"E:\Zuzka\CH_stacking");
	
	List<String> ligands = new List<String>();
	var files = File.ReadAllLines(@"LigandExpo_copy_structures.txt").ToHashSet(); //seznam ligandů ... už neaktuální (příště udělat přes složku a ne přes soubor)
	var locker = new object(); //potřeba pro paralelní běh
	var counter = 0; //počítadlo souborů, které už jsem prošla
	var ligandPaths = new List<String>(); //cesty k vybraným ligandům

	/*
    	definice kruhu pentosy a hexosy (pro předejítí nedorozuměním s Lukášem, pojmenovány tetrosa a pentosa)
		!!!! neobsahují podmínku o násobnosti vazby !!!!
	*/

	var tetroseC3 = QueryBuilder.AtomNames(new string[]{"C3"}).
						   Inside(QueryBuilder.Rings(new string[]{"C", "C", "C", "C", "O"})).
						   Filter(l => QueryBuilder.IsConnectedTo(l, QueryBuilder.Atoms(new string[]{"O"})));

	var tetroseC4 = QueryBuilder.AtomNames(new string[]{"C4"}).
						   Inside(QueryBuilder.Rings(new string[]{"C", "C", "C", "C", "O"})).
						   Filter(l => QueryBuilder.IsConnectedTo(l, QueryBuilder.Atoms(new string[]{"O"})));
						   

	var pentoseC3 = QueryBuilder.AtomNames(new string[]{"C3"}).
						   Inside(QueryBuilder.Rings(new string[]{"C", "C", "C", "C", "C", "O"})).
						   Filter(l => QueryBuilder.IsConnectedTo(l, QueryBuilder.Atoms(new string[]{"O"})));
						   
	var pentoseC4 = QueryBuilder.AtomNames(new string[]{"C4"}).
						   Inside(QueryBuilder.Rings(new string[]{"C", "C", "C", "C", "C", "O"})).
						   Filter(l => QueryBuilder.IsConnectedTo(l, QueryBuilder.Atoms(new string[]{"O"})));

	var tetrose = QueryBuilder.Or(tetroseC3, tetroseC4);
	var pentose = QueryBuilder.Or(pentoseC3, pentoseC4);
	var query = QueryBuilder.Or(tetrose, pentose).ToMetaQuery().Compile();	

	Parallel.ForEach(files, new ParallelOptions(){MaxDegreeOfParallelism = 4}, e =>
	{
		try{
			var str = StructureReader.Read(e).Structure;  // nacti strukturu
	
			var result = query.Matches(str); // podivej se, jestli ma query ve strukture nejaky match
									 //v result jsou uloženy všechny patterny pro jednu strukturu      
	
			if(result.Count != 0) {
				ligands.Add(str.Id); 
				ligandPaths.Add(e);
			}
	
			lock (locker){	
				counter++;
				counter.Dump(); //pro přehled během běhu programu
			}
	
		}catch(Exception ex){
			ex.Dump();
		}
	});

	File.WriteAllLines("ligandPathsC3C4.txt", ligandPaths);
	File.WriteAllLines("ligandsC3C4.txt", ligands);
}