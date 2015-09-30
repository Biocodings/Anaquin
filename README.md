Anaquin
=======

Anaquin is a C++/R bioinformatics framework for spiked-ins that provide qualitative and quantitative control for next-generation sequencing. 

The project was started in 2015 by <a href='https://www.google.com.au/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=0CB4QFjAAahUKEwiWt5b-7p3IAhWEjJQKHcxhDMg&url=http%3A%2F%2Fwww.garvan.org.au%2F&usg=AFQjCNF03pFvjJsIYqEbmxMV3SBTC5PJxg&sig2=jxHlEHfy_CNSJ4cZyVfvVQ'>Garvan Institute of Medical Research</a>. See <a href='http://www.anaquin.org'>www.anaquin.org</a> for further details. In particular, the <a href='http://www.anaquin.org/about/introduction/'>introduction page</a> should be helpful.

The project is maintained by <b>Ted Wong</b>, binformatics software engineer at Garvan Institute.

Overview
========

<img src='http://www.anaquin.org/wp-content/uploads/2015/07/About_Intro_F1_3.svg'>

Next-generation sequencing (NGS) enables rapid, cheap and high-throughput determination of DNA (or RNA) sequences within a user’s sample. NGS methods have been applied widely, and have fuelled major advances in the life sciences and clinical health care over the past decade. However, NGS typically generates a large amount of sequencing data that must be first analyzed and interpreted with bioinformatic tools. There is no standard way to perform an analysis of NGS data; different tools provide different advantages in different situations. For example, the tools for finding small deletions in the human genome are different to the tools to find large deletions. The sheer number, complexity and variation of sequences further compound this problem, and there is little reference by which compare next-generation sequencing and analysis.

To address this problem, we have developed a suite of synthetic nucleic-acid standards that we term sequins. Sequins represent genetic features, such as genes, large structural rearrangements, that are often analyzed with NGS. However, whilst sequins may act like a natural genetic feature, their primary sequence is artificial, with no extended homology to natural genetic sequences. Sequins are fractionally added to the extracted nucleic-acid sample prior to library preparation, so they are sequenced along with your sample of interest. The reads that derive from sequins can be identified by their artificial sequences that prevent their cross-alignment to the genome of known organisms. 

