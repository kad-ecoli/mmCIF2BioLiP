const navLinks = [
  { name: "D-I-TASSER", url: "/D-I-TASSER" },
  { name: "DMFold", url: "/DMFold" },
  { name: "HPmod", url: "/HPmod" },
  { name: "LOMETS3", url: "/LOMETS" },
  { name: "StarFunc", url: "/StarFunc/" },
  { name: "InterLabelGO+", url: "/InterLabelGO/" },
  { name: "US-align", url: "/US-align" },
  { name: "DeepMSA2", url: "/DeepMSA" },
  { name: "BioLiP2", url: "/BioLiP" },
  { name: "FURNA", url: "/furna" },
  { name: "ShapeME", url: "/shapeme" }
];

const navElement = document.getElementById('sidebar');
const list = document.createElement('ul');

navLinks.forEach(link => {
  const li = document.createElement('li');
  const a = document.createElement('a');
  a.href = link.url;
  a.textContent = link.name;
  li.appendChild(a);
  list.appendChild(li);
});

navElement.appendChild(list);
