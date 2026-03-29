let RDKit;

initRDKitModule().then(instance => {
    RDKit = instance;
});


document.querySelectorAll(".tab").forEach(tab => {
  tab.addEventListener("click", () => {
    const mode = tab.dataset.mode;

    document.querySelectorAll(".tab").forEach(t => t.classList.remove("active"));
    tab.classList.add("active");

    document.querySelectorAll(".tab-content").forEach(c => c.classList.remove("active"));
    document.getElementById("mode-" + mode).classList.add("active");

    if (mode === "advanced") {
      document.getElementById("layers-compare").style.display = "none";
      document.getElementById("layers-advanced").style.display = "block";
    } else {
      document.getElementById("layers-compare").style.display = "block";
      document.getElementById("layers-advanced").style.display = "none";
    }
  });
});

function updateLayerVisibility(mode) {
    const container = document.getElementById("layers-container");
    if (mode === "advanced") {
        container.querySelectorAll("input.level-checkbox").forEach(cb => cb.disabled = false);
    } else {
        container.querySelectorAll("input.level-checkbox").forEach(cb => cb.disabled = true);
    }
}

async function draw(inchi1, inchi2) {
    visualizeFromInchi("mol1", inchi1);
    visualizeFromInchi("mol2", inchi2);
}

async function visualizeFromInchi(containerId, inchi) {
    const element = document.getElementById(containerId);
    element.innerHTML = "";

    try {
        const res = await fetch("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/SDF?record_type=3d", {
            method: "POST",
            headers: { "Content-Type": "application/x-www-form-urlencoded" },
            body: `inchi=${encodeURIComponent(inchi)}`
        });

        if (!res.ok) throw new Error();
        const data = await res.text();

        const viewer = $3Dmol.createViewer(element, { backgroundColor: 'white' });
        viewer.addModel(data, "sdf");
        viewer.setStyle({}, { stick: { radius: 0.15 }, sphere: { scale: 0.25 } });
        viewer.zoomTo();
        viewer.render();

    } catch {
        try {
            const mol = RDKit.get_mol_from_inchi(inchi);
            element.innerHTML = mol.get_svg();
        } catch {
            element.innerHTML = "Invalid InChI";
        }
    }
}

async function compare(isAdvanced = false) {
    const inchi1 = isAdvanced ? document.getElementById("inchi1_adv").value.trim() : document.getElementById("inchi1").value.trim();
    const inchi2 = isAdvanced ? document.getElementById("inchi2_adv").value.trim() : document.getElementById("inchi2").value.trim();

    if (!inchi1 || !inchi2) {
        showToast("Please enter both InChIs", "error");
        return;
    }

    draw(inchi1, inchi2);

    let url = "http://127.0.0.1:5000/api/compare_inchis";
    let body = { inchi1, inchi2 };

    if (isAdvanced) {
        url = "http://127.0.0.1:5000/api/compare_inchis_custom";
        body.levels = Array.from(document.querySelectorAll(".level-checkbox:checked")).map(cb => cb.value);
        console.log("SELECTED LEVELS:", body.levels);
    }

    try {
        const res = await fetch(url, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(body)
        });

        const text = await res.text();

        console.log("RAW RESPONSE:", text);

        let data;
        try {
            data = JSON.parse(text);
        } catch {
            showToast("Server returned invalid JSON", "error");
            return;
        }

        if (!res.ok) {
            showToast(data.message || "Error", "error");
            return;
        }
        if (!res.ok) { showToast(data.message || "Error", "error"); return; }

        updateLayers(mapBackendResults(data.results), isAdvanced);

    } catch (err) {
        console.error(err);
        showToast("Server error", "error");
    }
}

function mapBackendResults(results) {
    return {
        complete_identity: results.COMPLETE_IDENTITY,
        isotope: results.ISOTOPIC_INDEPENDENCE,
        salt: results.SALTS_INDEPENDENCE,
        charge: results.CHARGES_INDEPENDENCE,
        stereo_cis_trans: results.STEREOCHEMICAL_CIS_TRANS_INDEPENDENCE,
        double_bond: results.DOUBLE_BONDS_INDEPENDENCE,
        tautomer: results.TAUTOMER_INDEPENDENCE,
        substituent: results.SUBSTITUENT_POSITION_INDEPENDENCE
    };
}

function updateLayers(results, isAdvanced = false) {
    document.querySelectorAll(".layer").forEach(layer => {

        let key = null;
        let checkbox = layer.querySelector("input");

        if (checkbox) {
            key = checkbox.value;

            if (isAdvanced && !checkbox.checked) {
                layer.style.display = "none";
                return;
            } else {
                layer.style.display = "block";
            }

        } else {
            key = layer.dataset.key;
        }

        if (!key) return;

        const match = results[key];

        let badge = layer.querySelector(".badge");
        if (!badge) {
            badge = document.createElement("span");
            badge.classList.add("badge");
            layer.appendChild(badge);
        }

        if (match === true) {
            layer.classList.add("match");
            layer.classList.remove("nomatch");
            badge.className = "badge green";
            badge.innerText = "INDEPENDENT";
        } else if (match === false) {
            layer.classList.add("nomatch");
            layer.classList.remove("match");
            badge.className = "badge red";
            badge.innerText = "NOT INDEPENDENT";
        } else {
            layer.classList.remove("match", "nomatch");
            badge.className = "badge";
            badge.innerText = "N/A";
        }
    });
}


function showToast(message, type = "info", duration = 4000) {
    const container = document.getElementById("toast-container");
    if (!container) return;

    const toast = document.createElement("div");
    toast.classList.add("toast", `toast-${type}`);
    toast.textContent = message;

    container.appendChild(toast);

    setTimeout(() => {
        toast.style.transition = "transform 0.3s ease, opacity 0.3s ease";
        toast.style.transform = "translateX(100%)";
        toast.style.opacity = "0";

        toast.addEventListener("transitionend", () => toast.remove());
    }, duration);
}