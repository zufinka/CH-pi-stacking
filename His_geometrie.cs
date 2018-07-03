using System;
using WebChemistry.Framework.Core;
using WebChemistry.Queries.Core;
using WebChemistry.Framework.Core.Pdb;
using WebChemistry.Framework.Math;

/*
    protřídění vybraných motivů (v .pdb formátu) podle geometrických parametrů
    v případě Trp zjišťovány dodatečné vlastnosti
 */
void Main()
{
	Directory.SetCurrentDirectory(@"E:\Zuzka\CH_stacking");//výběr pracovního adresáře
	
	var files = Directory.GetFiles("motivy/C3C4/HIS/all/vodiky/ligand"); //seznam PDB motivů s vodíky
    files.Count().Dump();//pro přehled
	
	/*
        patternQuery dotazy:
        vyfiltrování motivů se správnou vzdáleností a zároveň uložení atomů, které jsou potřeba pro torzní úhel (CH atomy aromatického kruhu + CH atomy napojené na kruh)
	*/								
	//HIS
	var stacking2 = QueryBuilder.Cluster(4.5, 
						QueryBuilder.RingAtoms(QueryBuilder.Atoms(new string[]{"C"})).Inside(QueryBuilder.Residues(new string[]{"HIS"})).Union(),
						QueryBuilder.Or(QueryBuilder.Rings(new string[]{"C", "C", "C", "C", "O"}), QueryBuilder.Rings(new string[]{"C", "C", "C", "C", "C", "O"})).ConnectedAtoms(1).
							Flatten(a => QueryBuilder.Find(a, QueryBuilder.Atoms(new string[]{"C"}).
								Filter(c => QueryBuilder.IsConnectedTo(c, QueryBuilder.Atoms(new string[]{"H"}))))).Union()).
									Flatten(l => QueryBuilder.Find(l, QueryBuilder.Atoms())).ToMetaQuery().Compile();
									
    /*
        proměnné potřebné pro zápis do souboru
     */
    var zapisSet = new HashSet<String>();
	zapisSet.Add("name;ligand;stacking distance;stacking torsion angle;signature;chain;");
    var ligandName = new List<String>();//ligandy ve vybraných motivech
    var pdbNames = new HashSet<String>();//PDB ID struktur, ve kterých byly nalezeny hledané motivy
    int motiveCount = 1;//počítadlo motivů v rámci jedné PDB struktury
    var lastPdbID = "";//proměnná pro ošetření situace u první struktury, která se prochází
    var motivesCount = new Dictionary<String, int>();//ukládá celkový počet motivů u jedné struktury (PDB ID)
    var pdbLigands = new Dictionary<String, List<String>>();//ukládá seznam ligandů k jednotlivým strukturám (PDB ID)
    var ligandsHistogram = new Dictionary<String, int>();//ukládá seznam ligandů a jejich počet
	
	foreach (var e in files)
	{
		var str = StructureReader.Read(e).Structure;
		
		//HIS řazení: CD2, CE1, CG, atomy cukru
        /*
            varianty proměnné sorted, podle toho, s čím aktuálně potřebuju pracovat
         */
		//var sorted = stacking2.Matches(str).OrderBy(x => x.Atoms.First().IsHetAtom()).ThenBy(x => x.Atoms.First().PdbName()).Select(x => x.Atoms.First().Position).ToList();
		//var sorted = stacking2.Matches(str).OrderBy(x => x.Atoms.First().IsHetAtom()).ThenBy(x => x.Atoms.First().PdbName()).ToList();
		var sorted = stacking2.Matches(str).OrderBy(x => x.Atoms.First().IsHetAtom()).ThenBy(x => x.Atoms.First().PdbName()).Select(x => x.Atoms).ToList();
		
        /*
            projdu případy, kdy se našly motivy s aromatickou amk v blízkosti cukru
            aromatická amk bude nalezená vždy!!!
         */
		if (sorted.Count > 3) //závisí na počtu atomů arom. residua, které ukládám
		{	
			//HIS: výpočet středu kruhu = (střed úsečky CD2 a CG) a CE1
			double cdg2X = (sorted.ElementAt(0).First().Position.X + sorted.ElementAt(2).First().Position.X) / 2;
			double cdg2Y = (sorted.ElementAt(0).First().Position.Y + sorted.ElementAt(2).First().Position.Y) / 2;
			double cdg2Z = (sorted.ElementAt(0).First().Position.Z + sorted.ElementAt(2).First().Position.Z) / 2;
			double aromCenterX = (cdg2X + sorted.ElementAt(1).First().Position.X) / 2;
			double aromCenterY = (cdg2Y + sorted.ElementAt(1).First().Position.Y) / 2;
			double aromCenterZ = (cdg2Z + sorted.ElementAt(1).First().Position.Z) / 2;
			Vector3D aromCenter = new Vector3D(aromCenterX, aromCenterY, aromCenterZ);

            /*
                výpočet vzdálenosti mezi nejbližším CH atomem cukru a středem aromatické části
             */
			var nejblizsiSouradnice = new Vector3D();
			String nejblizsiResiduum = "";
			double minDistance = 10;

			for (int i = 3; i < sorted.Count; i++)//hledám nejbližší atom cukru = procházím až atomy cukru !!!zkontrolovat kolik atomů aromatické amk je před tím!!!
			{	
				//TYR + PHE + HIS + TRP na části
				Vector3D d = new Vector3D( (sorted.ElementAt(i).First().Position.X - aromCenter.X), 
											(sorted.ElementAt(i).First().Position.Y - aromCenter.Y), 
											(sorted.ElementAt(i).First().Position.Z - aromCenter.Z) );
				
				var distance = d.Length;
				
				if ((distance < minDistance) && (sorted.ElementAt(i).First().PdbResidueName() != "HIS") )//teoreticky můžu mít v sorted víc aromatických residuí než jedno
				{
					nejblizsiSouradnice = sorted.ElementAt(i).First().Position;
					nejblizsiResiduum = sorted.ElementAt(i).First().PdbResidueName();
					minDistance = distance;
				}
			}
			
			/*
                výpočet úhlu pro filtr
            */
            //rovina aromatického kruhu
			Plane3D rovina = Plane3D.FromPoints(sorted.ElementAt(0).First().Position, sorted.ElementAt(1).First().Position, aromCenter);
			
            //vzdálenost nejbližšího CH atomu od roviny aromatického kruhu
            var cPlaneDistance = Math.Abs( rovina.A*nejblizsiSouradnice.X + rovina.B*nejblizsiSouradnice.Y + rovina.C*nejblizsiSouradnice.Z + rovina.D )/
                                    Math.Sqrt( Math.Pow(rovina.A, 2) + Math.Pow(rovina.B, 2) + Math.Pow(rovina.C, 2) );

			double cosTor1 = cPlaneDistance / minDistance;
								 
			var filterAngleRad = Math.Acos(cosTor1);
			var filterAngle = filterAngleRad * (180 / Math.PI);
			
			//uložit: označení motivu; vzdálenost cukr-arom. amk; torzni uhel
			if ((minDistance < 4.5) && (filterAngle < 45)) //filtr pro výběr motivů
			{
				var zapis = new StringBuilder();
				zapis.Append(Path.GetFileNameWithoutExtension(e));
				zapis.Append(";");
				zapis.Append(nejblizsiResiduum);
				zapis.Append(";");
				zapis.Append(minDistance);
				zapis.Append(";");
				zapis.Append(filterAngle);
				zapis.Append(";");
				foreach (var res in str.PdbResidues()){
					zapis.Append(res.Name);
					zapis.Append(" ");
				}
				zapis.Append(";");
				foreach (var chain in str.PdbChains()){
					zapis.Append(chain.Key);
				}
				
				zapis.Append(";");
                ligandName.Add(nejblizsiResiduum);

                //zápis vybraného motivu do souboru
				var path = new StringBuilder("motivy/C3C4/HIS/5_11_2017/ligand/");
				path.Append(Path.GetFileNameWithoutExtension(e));
				using (TextWriter write = File.CreateText(path.ToString())){
					str.WritePdb(write);
				}

				zapisSet.Add(zapis.ToString());
				
				var pdbID = Path.GetFileNameWithoutExtension(e).Substring(0,4);
				if (pdbNames.Add(pdbID) == false)
				{
					motiveCount += 1;
					pdbLigands[pdbID].Add(nejblizsiResiduum);
				}else
				{
					if (lastPdbID != "")
					{
						motivesCount.Add(lastPdbID, motiveCount);
					}
					pdbLigands.Add(pdbID, new List<String>());
					pdbLigands[pdbID].Add(nejblizsiResiduum);
					motiveCount = 1;
					lastPdbID = pdbID;
				}
			}
		}
	}
	
	motivesCount.Add(lastPdbID, motiveCount);

    File.WriteAllLines("HIS_ligand.csv", zapisSet);
	
    //výpočet histogramu jednotlivých ligandů
    var ligandNames = new HashSet<String>();//pomocná proměnná pro zjištění, jestli už jsem daný ligand měla nebo ne
	foreach (var lig in ligandName)
	{
		if (ligandNames.Add(lig) == false)
		{
			ligandsHistogram[lig] += 1;
		}else
		{
			ligandsHistogram.Add(lig, 1);
		}
    }
	
	int counter = 1;
	using (TextWriter writer = File.CreateText("HIS_ligand_pdb.csv"))
	{
		writer.WriteLine("pdbID;motives count;ligands");
		foreach (var a in pdbNames)
		{
			counter = 1;
			writer.Write(a);
			writer.Write(";");
			writer.Write(motivesCount[a]);
			writer.Write(";");
			foreach (var b in pdbLigands[a])
			{
				if (counter != 1)
				{
					writer.Write(",");
				}
				writer.Write(b);
				counter += 1;
			}
			writer.Write(";");
			writer.Write(writer.NewLine);
		}
	}
	
	using (TextWriter wr = File.CreateText("HIS_ligand_ligands.csv"))
	{
		wr.WriteLine("ligand;count;");
		foreach (var l in ligandNames)
		{
			wr.Write(l);
			wr.Write(";");
			wr.Write(ligandsHistogram[l]);
			wr.Write(";");
			wr.Write(wr.NewLine);
		}
	}
}