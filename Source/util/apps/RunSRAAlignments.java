package util.apps;
import java.io.*;
import util.gen.*;

public class RunSRAAlignments {


	public static void main(String[] args) {
		if (args.length == 0) Misc.printExit("Enter the directory containing xxx.lite.sra files and a directory in which to save the scripts.");
		
		//pull files
		File dir = new File(args[0]);
		File[] sras = IO.fetchFilesRecursively(dir, ".lite.sra");
		File scripts = new File(args[1]);
		
		//make and save scripts
		for (int i=0; i< sras.length; i++){
			//extract name
			String name = sras[i].getName().replace(".lite.sra", "");
			//fetch script
			String script = fetchScript(name);
			//save
			IO.writeString(script, new File(scripts, name+".sh"));
		}

	}
	
	public static String fetchScript(String name){
		return
		"#PBS -l walltime=24:00:00,nodes=1:ppn=12 \n"+
		"#PBS -m a \n"+
		"#PBS -M david.nix@hci.utah.edu \n"+
		"#PBS -N "+name+" \n"+
		"#PBS -j oe \n"+
		"#PBS -A kaplan \n"+
		"#PBS -o /uufs/chpc.utah.edu/common/home/u0028003/Scratch/Run3/Logs/"+name+".log \n"+
		" \n"+
		"#set params \n"+
		"n="+name+" \n"+
		" \n"+
		"s=/scratch/local/u0028003/SRAs \n"+
		"f=/scratch/local/u0028003/Fastq \n"+
		"a=/scratch/local/u0028003/Alignments \n"+
		"novoalign=/uufs/chpc.utah.edu/common/home/u0028003/BioApps/Novocraft/novocraft/novoalign \n"+
		"sra=/uufs/chpc.utah.edu/common/home/u0028003/BioApps/sratoolkit.2.0b5-centos_linux64 \n"+
		"param='-rRandom -h100 -t60 -b2' \n"+
		"index=/scratch/local/u0028003/hg18AdapterPhixLambdaBisulfite.novoindex \n"+
		"homeIndex=/uufs/chpc.utah.edu/common/home/u0028003/Genomes/hg18AdapterPhixLambdaBisulfite.novoindex \n"+
		"homeSRAs=/uufs/chpc.utah.edu/common/home/u0028003/Scratch/Run3/SRAs \n"+
		"homeAlignments=/uufs/chpc.utah.edu/common/home/u0028003/Scratch/Run3/Alignments \n"+
		"u=u0028003 \n"+
		" \n"+
		"echo Print name and date \n"+
		"echo $n \n"+
		"date \n"+
		" \n"+
		"echo Change default file permissions \n"+
		"umask 000 \n"+
		" \n"+
		"echo Delete everything in temp \n"+
		"rm -rf /scratch/local/* \n"+
		" \n"+
		"echo Make tempDirs \n"+
		"mkdir /scratch/local/$u \n"+
		"mkdir $s \n"+
		"mkdir $f \n"+
		"mkdir $a \n"+
		" \n"+
		"echo Copy over files from ibrix \n"+
		"cp $homeSRAs/$n'.lite.sra' $s/ \n"+
		"cp $homeIndex $index \n"+
		" \n"+
		"echo Extract SRA data \n"+
		"$sra/fastq-dump -A $n -D $s/$n'.lite.sra' -O $f -F \n"+
		"echo \n"+
		" \n"+
		"echo Novoalignments \n"+
		"f0=$f/$n'.fastq' \n"+
		"f1=$f/$n'_1.fastq' \n"+
		"f2=$f/$n'_2.fastq' \n"+
		"rs=$a/$n'_Single.novo' \n"+
		"rp=$a/$n'_Pair.novo' \n"+
		"echo Number of lines in fastq \n"+
		"wc -l $f1 \n"+
		"$novoalign -d $index $param -f $f0 | grep chr > $rs \n"+
		"$novoalign -d $index $param -f $f1 $f2 | grep chr > $rp \n"+
		"echo \n"+
		" \n"+
		"echo Compress alignments \n"+
		"gzip $rs \n"+
		"gzip $rp \n"+
		"echo \n"+
		" \n"+
		"echo Move results to ibrix \n"+
		"mv $a/* $homeAlignments \n"+
		" \n"+
		"echo Clean up \n"+
		"rm -rf /scratch/local/$u \n"+
		" \n"+
		"echo \n"+
		"date \n"+
		"echo Done! \n";

	}

}
