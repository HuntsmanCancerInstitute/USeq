package util.bio.annotation;

import java.util.Comparator;

public class ComparatorConcatinatedSequence implements Comparator<ConcatinatedSequence> {

		
			/**Sorts by sequence*/
			public int compare(ConcatinatedSequence first, ConcatinatedSequence second) {
				String seqFirst = first.getSequence();
				String seqSecond = second.getSequence();
				return seqFirst.compareTo(seqSecond);
			}

		}


