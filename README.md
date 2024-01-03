
<h1 align="center">rNMP Hotpsots</h1>
To obtain rNMP hotspots(locations of highly abundant rNMPs) using a percentage or poisson pvalue threshold and map on respective reference genomes. 

<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--
-->
[![Commits][Commits-shield]][Commits-url]
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Website][website-shield]][website-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="##Installation">Installation</a></li>
      <ul>
        <li><a href="###Getting-the-code">Getting the code</a></li>
        <li><a href="###Creating-the-enviroment-with-required-dependencies">Creating the enviroment with required dependencies</a></li>
        <li><a href="###Additional-Dependencies">Additional Dependencies</a></li>
      </ul>
    </li>
    <li><a href="##Usage">Usage</a></li>
      <ul>
        <li><a href="###Defining-variables">Defining variables</a></li>
        <li><a href="###Normalization of bed files for coverage">Normalization of bed files for coverage</a></li>
        <li><a href="###Getting-hotspots-using-threshold">Getting counts using bed file</a></li>
       <li><a href="###Visualization">Visualization</a></li>
      </ul>
    <li><a href="##Contributing">Contributing</a></li>
    <li><a href="##License">License</a></li>
    <li><a href="##Contact">Contact</a></li>
    <li><a href="##Citations">Citations</a></li>
  </ol>
</details>

<!-- Installation -->
## Installation
### Getting the code
The development version from [GitHub](https://github.com/) with:
```sh
git clone https://github.com/DKundnani/rNMP_hotspots.git
```
### Creating the enviroment with required dependencies
```sh
conda env create --name RibosemapQC_env --file /rNMP_hotspots/yml/r_env.yml
```
### Additional Dependencies
* Reference genome files (.fa and .fai) of the organism being used(Also used to generate bed files)

<!-- USAGE -->

## Usage
="###Defining-variables
```bash
lib=path/to/AGS/ribo-DNA-order #First col FScode, 3rd col basename of bam files
bed=path/to/AGS/bed
dna=path/to/AGS/DNAseq/aligned
normbed=path/to/AGS/norm_counts
script=path/to/AGS/rNMP_hotspots
genome=path/to/reference/sacCer2/sacCer2-nucl.fa.fai![image](https://github.com/DKundnani/rNMP_hotspots/assets/20824535/9cbc97dc-9292-47ad-805d-e52b633e7447)
```
### Normalization of bed files for coverage
```bash

```
### Getting-hotspots-using-threshold
```bash

```
### Getting-hotspots-using-threshold
```bash

```
### Getting-hotspots-using-threshold
```bash

```


<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact
Deepali L. Kundnani - [deepali.kundnani@gmail.com](mailto::deepali.kundnani@gmail.com)    [![LinkedIn][linkedin-shield]][linkedin-url] 
<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ACKNOWLEDGMENTS -->
## Citations
Use this space to list resources you find helpful and would like to give credit to. I've included a few of my favorites to kick things off!+
* <b>Light-strand bias and enriched zones of embedded ribonucleotides are associated with DNA replication and transcription in the human-mitochondrial genome. </b>
Penghao Xu, Taehwan Yang, Deepali L Kundnani, Mo Sun, Stefania Marsili, Alli L Gombolay, Youngkyu Jeon, Gary Newnam, Sathya Balachander, Veronica Bazzani, Umberto Baccarani, Vivian S Park, Sijia Tao, Adriana Lori, Raymond F Schinazi, Baek Kim, Zachary F Pursell, Gianluca Tell, Carlo Vascotto, Francesca Storici
<i>  Acids Research </i> 2023;, gkad1204, [https://doi.org/10.1093/nar/gkad1204](https://doi.org/10.1093/nar/gkad1204)
* <b>Distinct features of ribonucleotides within genomic DNA in Aicardi-Goutières syndrome (AGS)-ortholog mutants of <i>Saccharomyces cerevisiae</i> </b>
Deepali L. Kundnani, Taehwan Yang, Alli L. Gombolay, Kuntal Mukherjee, Gary Newnam, Chance Meers, Zeel H. Mehta, Celine Mouawad, Francesca Storici
bioRxiv 2023.10.02.560505; doi:[https://doi.org/10.1101/2023.10.02.560505]( https://doi.org/10.1101/2023.10.02.560505)
* Kundnani, D. (2024). Ribose-Map-QC (1.0.0). Zenodo. [https://doi.org/10.5281/zenodo.10455801](https://doi.org/10.5281/zenodo.10455801) [![DOI](https://zenodo.org/badge/382943161.svg)](https://zenodo.org/doi/10.5281/zenodo.10453008)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/DKundnani/rNMP_hotspots?style=for-the-badge
[contributors-url]: https://github.com/DKundnani/rNMP_hotspots/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/DKundnani/rNMP_hotspots?style=for-the-badge
[forks-url]: https://github.com/DKundnani/rNMP_hotspots/forks
[stars-shield]: https://img.shields.io/github/stars/DKundnani/rNMP_hotspots?style=for-the-badge
[stars-url]: https://github.com/DKundnani/rNMP_hotspots/stargazers
[issues-shield]: https://img.shields.io/github/issues/DKundnani/rNMP_hotspots?style=for-the-badge
[issues-url]: https://github.com/DKundnani/rNMP_hotspots/issues
[license-shield]: https://img.shields.io/github/license/DKundnani/rNMP_hotspots?style=for-the-badge
[license-url]: https://github.com/DKundnani/rNMP_hotspots/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/deepalik
[product-screenshot]: images/screenshot.png
[commits-url]: https://github.com/DKundnani/rNMP_hotspots/pulse
[commits-shield]: https://img.shields.io/github/commit-activity/t/DKundnani/rNMP_hotspots?style=for-the-badge
[website-shield]: https://img.shields.io/website?url=http%3A%2F%2Fdkundnani.bio%2F&style=for-the-badge
[website-url]:http://dkundnani.bio/ 
