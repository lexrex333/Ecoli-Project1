# Ecoli-Project1
## Design Doc
### Overview
Escherichia coli is a major component that surrounds life on Earth, varying from causing disease to providing benefits within human microbiomes. In just the human colon alone, there are around a billion E. coli cells (LTEE Introduction)! Research of E. coli has allowed for significant discoveries in the dynamics and attributes of microbial life. E. coli plays a pivotal role in this research as it can proliferate rapidly. The bacteria only take a few hours to grow 100-fold from its initial population size, representing 6 and ⅔ generations. Because of this rapid duplication that occurs in bacteria, horizontal gene transfer and genome rearrangement occur at a significantly greater rate than in other types of populations. The variation between different generations of bacteria can develop quickly as the proliferation time is minimal. As a result, the documentation of each generation can provide insight into the development of variation. The Long Term Evolution Experiment is a great example of the documentation of E. coli, especially since E. coli can be preserved and used at a later date.

The Long Term Evolution Experiment observes and quantifies the process of evolution in action by growing twelve populations of E. coli that were subcultured from two almost-identical ancestral strains. These two strains, REL606 and REL607, are each ancestors to six ongoing populations. Every day, each population was propagated by transferring 1% of the previous day’s culture into fresh medium. The glucose medium acts as a limiting nutrient for the strains, and runs out once they have grown 100-fold, and propagation is repeated the next day. 

This project aims to compare the evolutionary divergence over time within each of these twelve E. coli populations throughout 50,000 generations. Using 264 genomes sequenced by Tenaillon et al., we compare each population to its ancestral strain (REL606 or REL607) from 1988. Here, we present a pipeline to quantify and visualize the genome similarities of each strain in relation to its ancestors using metrics of average nucleotide identity (ANI), Mash, and digital DNA-DNA hybridization (dDDH). Previous research from Tenaillon et al. observed genetic changes within each population over time, but did not use these same metrics (ANI, Mash, dDDH) to compare populations to their ancestors.

##### References: 
[LTEE Introduction](https://the-ltee.org/about/)

[Tenaillon et al.](https://www.nature.com/articles/nature18959)
