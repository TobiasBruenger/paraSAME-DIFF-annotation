<h1>Annotating para-PS1/PM5 Criteria with a Custom Variant Dataset for Gene Families</h1>

<p>Welcome to the GitHub repository associated with the manuscript titled "Conserved Missense Variant Pathogenicity and Correlated Phenotypes across Paralogous Genes," which is currently under submission. Upon publication, we will provide a reference link to the article.</p>

<h2>Overview</h2>

<p>This repository contains R scripts that empower users to perform the following tasks:</p>

<h3>1. Generate a Variant Evidence Table</h3>

<p>These R scripts enable users to create a table that identifies protein residues in genes from a selected gene family where the presence of a pathogenic variant in a paralogous gene may offer evidence for evaluating a novel variant found in an affected individual. To generate this table, users will need:</p>

<ul>
  <li>An alignment file for their gene family of interest. Pre-annotated alignment files can be obtained from <a href="https://github.com/dlal-group/PERs">this link</a>. Please acknowledge the original publication by Pérez-Palma et al. (2020) if you use these alignment files (<a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6961572/">link to the publication</a>).</li>
  <li>A custom set of pathogenic variants that will serve as evidence if they are present in a conserved manner and at the same protein alignment position in a paralogous gene.</li>
</ul>

<h3>2. Variant Pathogenicity Annotation</h3>

<p>The second script provided here takes the output from the first step or alternatively, you can use Supplementary Table 1 from the manuscript. This script annotates whether a custom set of variants provides evidence for the pathogenicity of a variant based on paralogous pathogenic variants.</p>

