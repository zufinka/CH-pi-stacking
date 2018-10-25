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
	var vodikova_vazba = QueryBuilder.Cluster(5, 
						QueryBuilder.Or(QueryBuilder.Atoms(new string[]{"O"}), QueryBuilder.Atoms(new string[]{"N"})).Inside(QueryBuilder.Residues(new string[]{"HIS"})).Union(),
						QueryBuilder.Or(QueryBuilder.Rings(new string[]{"C", "C", "C", "C", "O"}), QueryBuilder.Rings(new string[]{"C", "C", "C", "C", "C", "O"})).ConnectedAtoms(2).
							Flatten(a => QueryBuilder.Find(a, QueryBuilder.Or(QueryBuilder.Atoms(new string[]{"O"}), QueryBuilder.Atoms(new string[]{"N"})))).Union()).
									Flatten(l => QueryBuilder.Find(l, QueryBuilder.Atoms())).ToMetaQuery().Compile();
									
    /*
        proměnné potřebné pro zápis do souboru
     */
    var zapisSet = new HashSet<String>();
	zapisSet.Add("name;ligand;aminoacid_atom;sugar_atom;distance;angle;");
    var ligandName = new List<String>();//ligandy ve vybraných motivech
    var pdbNames = new HashSet<String>();//PDB ID struktur, ve kterých byly nalezeny hledané motivy
    int motiveCount = 1;//počítadlo motivů v rámci jedné PDB struktury
    var lastPdbID = "";//proměnná pro ošetření situace u první struktury, která se prochází
    var motivesCount = new Dictionary<String, int>();//ukládá celkový počet motivů u jedné struktury (PDB ID)
    var pdbLigands = new Dictionary<String, List<String>>();//ukládá seznam ligandů k jednotlivým strukturám (PDB ID)
    var ligandsHistogram = new Dictionary<String, int>();//ukládá seznam ligandů a jejich počet

	var pocetNepocitanychMotivu = 0; //pro kontrolu počtu motivů s chybou při výpočtu
	var pocet = 0; //počet nalezených motivů
	
	foreach (var e in files)
	{
		var str = StructureReader.Read(e).Structure;
		
		//HIS řazení: N, ND1, NE2, O !! pozor na řazení více HIS (N, N, ND1, ND1, NE2, NE2, O, O)
        /*
            varianty proměnné sorted, podle toho, s čím aktuálně potřebuju pracovat
         */
		//var sorted = vodikova_vazba.Matches(str).OrderBy(x => x.Atoms.First().IsHetAtom()).ThenBy(x => x.Atoms.First().PdbName()).Select(x => x.Atoms.First().Position).ToList();
		//var sorted = vodikova_vazba.Matches(str).OrderBy(x => x.Atoms.First().IsHetAtom()).ThenBy(x => x.Atoms.First().PdbName()).ToList();
		var sorted = vodikova_vazba.Matches(str).OrderBy(x => x.Atoms.First().IsHetAtom()).ThenBy(x => x.Atoms.First().PdbName()).Select(x => x.Atoms).ToList();
        //sorted.Dump();
		
        /*
            projdu případy, kdy se našly motivy s aromatickou amk v blízkosti cukru
            aromatická amk bude nalezená vždy!!! (možná už bezpředmětná podmínka, z historických důvodů ponechána)
         */
		var pocetAtomuJedneAMK = 4;
		var pocetAtomuAMK = pocetAtomuJedneAMK; //zkontrolovat, jestli jsem ji začala používat všude, kde to číslo používám
		if (sorted.Count > pocetAtomuAMK) //závisí na počtu atomů arom. residua, které ukládám
		{	
			//zjištění kolik mám AMK (hrubou silou, předpokládám, že nemám víc než 4)
			if (sorted.ElementAt(pocetAtomuAMK).First().PdbResidueName().Equals("HIS"))
			{
				if (sorted.ElementAt(pocetAtomuAMK*2).First().PdbResidueName().Equals("HIS"))
				{
					if (sorted.ElementAt(pocetAtomuAMK*3).First().PdbResidueName().Equals("HIS"))
					{
						pocetAtomuAMK = pocetAtomuAMK * 4;
					}else
					{
						pocetAtomuAMK = pocetAtomuAMK * 3;
					}
				}else
				{
					pocetAtomuAMK = pocetAtomuAMK * 2;
				}
			}

			/*
                výpočet vzdálenosti mezi nejbližším atomem aminokyseliny a atomem cukru
             */
			var nejblizsiSouradnice = new Vector3D();
			String nejblizsiResiduum = "";
			double minDistance = 10;
            var nejblizsiAtomCukr = sorted.ElementAt(pocetAtomuAMK).First();
            var nejblizsiAtomAmk = sorted.ElementAt(0).First();
            var angle = -10.0;

            //porovnám každý N nebo O z aminokyseliny s každým N nebo O cukru
            for (int amk = 0; amk < sorted.Count; amk++)
            {
				if (sorted.ElementAt(amk).First().PdbResidueName().Equals("HIS"))
				{
                	for (int i = pocetAtomuAMK-1; i < sorted.Count; i++)//hledám nejbližší atom cukru = procházím až atomy cukru 
			   		{	
                    	if ((sorted.ElementAt(i).First().PdbResidueName() != "HIS") && !(sorted.ElementAt(i).First().ElementSymbol.ToString().Equals("C")))//atomy HIS už by tam být neměly, ale jsou tam...
						{
				    		Vector3D d = new Vector3D( (sorted.ElementAt(i).First().Position.X - sorted.ElementAt(amk).First().Position.X), 
													(sorted.ElementAt(i).First().Position.Y - sorted.ElementAt(amk).First().Position.Y), 
													(sorted.ElementAt(i).First().Position.Z - sorted.ElementAt(amk).First().Position.Z) );
				
				    		var distance = d.Length;
				
				    		if (distance < minDistance) //&& (sorted.ElementAt(amk).First().PdbResidueName().Equals("HIS")))
				    		{
					    		nejblizsiSouradnice = sorted.ElementAt(i).First().Position;
					    		nejblizsiResiduum = sorted.ElementAt(i).First().PdbResidueName();
					    		minDistance = distance;
                        		nejblizsiAtomCukr = sorted.ElementAt(i).First();
                        		nejblizsiAtomAmk = sorted.ElementAt(amk).First();
				    		}
						}
			    	}
				}
            }

            /*
                výpočet úhlu pro filtr
            */
			
            if (minDistance == 0.0)
			{
				"CHYBA".Dump(); //jen pro přehled, jak moc je problém s řazením atomů závažný
			}

            if (minDistance < 3.43 && minDistance != 0.0) //problém s řazením - atomy cukru se dostaly mezi atomy AMK (u HIS jen u jednoho motivu)
            {
                var angleQueryAMK = QueryBuilder.AtomIds(nejblizsiAtomAmk.Id).ConnectedAtoms(1).ToMetaQuery().Compile();
                var angleMotivAMK = angleQueryAMK.Matches(str).OrderBy(x => x.Atoms.First().IsHetAtom()).ThenBy(x => x.Atoms.First().Id).Select(x => x.Atoms).ToList();

				var angleQueryCukr = QueryBuilder.AtomIds(nejblizsiAtomCukr.Id).ConnectedAtoms(1).ToMetaQuery().Compile();
                var angleMotivCukr = angleQueryCukr.Matches(str).OrderBy(x => x.Atoms.First().IsHetAtom()).ThenBy(x => x.Atoms.First().Id).Select(x => x.Atoms).ToList();

				var vodik = false;
				//zjištění, jestli Náš atom na cukru má navázaný vodík
				for (var a = 0; a < angleMotivCukr.First().Count; a++)
				{
					if (angleMotivCukr.First().ElementAt(a).ElementSymbol.ToString().Equals("H"))
					{
						vodik = true;
					}
				}

				//atom, který má vodík musí být uprostřed při výpočtu úhlu
				if (vodik)
				{//uprostřed je O nebo N cukru
					for (var a = 0; a < angleMotivCukr.ElementAt(0).Count; a++)
					{
                    	//první pokus o výpočet
                        	/*var d1 = new Vector3D( (angleMotivCukr.ElementAt(0).ElementAt(a).Position.X - nejblizsiSouradnice.X),
                                                (angleMotivCukr.ElementAt(0).ElementAt(a).Position.Y - nejblizsiSouradnice.Y),
                                                (angleMotivCukr.ElementAt(0).ElementAt(a).Position.Z - nejblizsiSouradnice.Z) );
                        	var sin1 = Math.Abs(angleMotivCukr.ElementAt(0).ElementAt(a).Position.X - nejblizsiSouradnice.X) / d1.Length;

                        	var sin2 = Math.Abs(nejblizsiSouradnice.X - nejblizsiAtomAmk.Position.X) / minDistance;
                        	angle = (Math.Asin(sin1) + Math.Asin(sin2)) * (180 / Math.PI);*/
							
							//pokus o výpočet přes kosinovou větu
							var aa = new Vector3D( (angleMotivCukr.ElementAt(0).ElementAt(a).Position.X - nejblizsiAtomAmk.Position.X),
												   (angleMotivCukr.ElementAt(0).ElementAt(a).Position.Y - nejblizsiAtomAmk.Position.Y),
												   (angleMotivCukr.ElementAt(0).ElementAt(a).Position.Z - nejblizsiAtomAmk.Position.Z) );
												  
							var b = new Vector3D( (angleMotivCukr.ElementAt(0).ElementAt(a).Position.X - nejblizsiAtomCukr.Position.X),
												  (angleMotivCukr.ElementAt(0).ElementAt(a).Position.Y - nejblizsiAtomCukr.Position.Y),
												  (angleMotivCukr.ElementAt(0).ElementAt(a).Position.Z - nejblizsiAtomCukr.Position.Z) );
												  
							var c = new Vector3D( (nejblizsiAtomCukr.Position.X - nejblizsiAtomAmk.Position.X),
												  (nejblizsiAtomCukr.Position.Y - nejblizsiAtomAmk.Position.Y),
												  (nejblizsiAtomCukr.Position.Z - nejblizsiAtomAmk.Position.Z) );
												  
							var cos = (Math.Pow(2, b.Length) + Math.Pow(2, c.Length) - Math.Pow(2, aa.Length)) / (2 * b.Length * c.Length);
							angle = Math.Acos(cos) * (180 / Math.PI);

							if ((cos > 1) || (cos < -1))
							{
								"CHYBA3".Dump();
								pocetNepocitanychMotivu++;
							}
                	}
				}else
				{//uprostřed je N na AMK
					for (var a = 0; a < angleMotivAMK.ElementAt(0).Count; a++)
					{
                    	/*var d2 = new Vector3D( (angleMotivAMK.ElementAt(0).ElementAt(a).Position.X - nejblizsiAtomAmk.Position.X),
                                                (angleMotivAMK.ElementAt(0).ElementAt(a).Position.Y - nejblizsiAtomAmk.Position.Y),
                                                (angleMotivAMK.ElementAt(0).ElementAt(a).Position.Z - nejblizsiAtomAmk.Position.Z) );
                        	var sin1 = Math.Abs(angleMotivAMK.ElementAt(0).ElementAt(a).Position.X - nejblizsiAtomAmk.Position.X) / d2.Length;

                        	var sin2 = Math.Abs(nejblizsiSouradnice.X - nejblizsiAtomAmk.Position.X) / minDistance;
                        	angle = (Math.Asin(sin1) + Math.Asin(sin2)) * (180 / Math.PI);*/
							
							//pokus o výpočet přes kosinovou větu
							var b = new Vector3D( (angleMotivAMK.ElementAt(0).ElementAt(a).Position.X - nejblizsiAtomAmk.Position.X),
												  (angleMotivAMK.ElementAt(0).ElementAt(a).Position.Y - nejblizsiAtomAmk.Position.Y),
												  (angleMotivAMK.ElementAt(0).ElementAt(a).Position.Z - nejblizsiAtomAmk.Position.Z) );
												  
							var aa = new Vector3D( (angleMotivAMK.ElementAt(0).ElementAt(a).Position.X - nejblizsiAtomCukr.Position.X),
												  (angleMotivAMK.ElementAt(0).ElementAt(a).Position.Y - nejblizsiAtomCukr.Position.Y),
												  (angleMotivAMK.ElementAt(0).ElementAt(a).Position.Z - nejblizsiAtomCukr.Position.Z) );
												  
							var c = new Vector3D( (nejblizsiAtomCukr.Position.X - nejblizsiAtomAmk.Position.X),
												  (nejblizsiAtomCukr.Position.Y - nejblizsiAtomAmk.Position.Y),
												  (nejblizsiAtomCukr.Position.Z - nejblizsiAtomAmk.Position.Z) );
												  
							var cos = (Math.Pow(2, b.Length) + Math.Pow(2, c.Length) - Math.Pow(2, aa.Length)) / (2 * b.Length * c.Length);
							angle =  (Math.Acos(cos)) * (180 / Math.PI);

							if ((cos > 1) || (cos < -1))
							{
								"CHYBA3".Dump();
								pocetNepocitanychMotivu++;
							}
                	}
            }
			
			angle.Dump();
			if (angle == -10.0 && minDistance < 3.43)//z důvodů neošetření pár chyb - pro přehled, o kolik motivů jsem přišla
			{
				pocetNepocitanychMotivu++;
			}
			
			//uložit: označení motivu; vzdálenost cukr-arom. amk; uhel
			if ((angle > 89) && (angle < 181)) //filtr pro výběr motivů zatím vypnutý, jen ošetřuji případy, kdy se úhel nepočítal, protože minimální vzdálenost byla příliš veliká
			{
				pocet++;
				var zapis = new StringBuilder();
				zapis.Append(Path.GetFileNameWithoutExtension(e));
				zapis.Append(";");
				zapis.Append(nejblizsiResiduum);
				zapis.Append(";");
                zapis.Append(nejblizsiAtomAmk.ElementSymbol);
                zapis.Append(";");
				zapis.Append(nejblizsiAtomCukr.ElementSymbol);
				zapis.Append(";");
				zapis.Append(minDistance);
				zapis.Append(";");
				zapis.Append(angle);
				zapis.Append(";");
				/*foreach (var res in str.PdbResidues()){
					zapis.Append(res.Name);
					zapis.Append(" ");
				}
				zapis.Append(";");
				foreach (var chain in str.PdbChains()){
					zapis.Append(chain.Key);
				}
				
				zapis.Append(";");*/
                ligandName.Add(nejblizsiResiduum);

                //zápis vybraného motivu do souboru
				var path = new StringBuilder("motivy/C3C4/HIS/H_24_10_2018/ligand/");
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

	"chyba".Dump();
	pocetNepocitanychMotivu.Dump();
	"nalezeno".Dump();
	pocet.Dump();

	motivesCount.Add(lastPdbID, motiveCount);

    File.WriteAllLines("HIS_H_ligand.csv", zapisSet);
	
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
	using (TextWriter writer = File.CreateText("HIS_H_ligand_pdb.csv"))
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
	
	using (TextWriter wr = File.CreateText("HIS_H_ligand_ligands.csv"))
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