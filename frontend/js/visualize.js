var RDKit = null;
var _rdkitReady = initRDKitModule().then(function(inst) { RDKit = inst; }).catch(function() {});
var _sdfCache = {};
var _fetching = {};

var _modalBackdrop   = null;
var _modalInchiLabel = null;
var _modalCanvas     = null;
var _modalViewer     = null;
var _modalBuilt      = false;

function visualizeFromInchi(containerId, inchi) {
    var el = document.getElementById(containerId);
    if (!el || !inchi) return;
    _renderCard(el, inchi);
}

function draw(inchi1, inchi2, id1, id2) {
    visualizeFromInchi(id1 || "mol1", inchi1);
    visualizeFromInchi(id2 || "mol2", inchi2);
}

function drawPair(leftEl, rightEl, inchi1, inchi2) {
    if (leftEl  && inchi1) _renderCard(leftEl,  inchi1);
    if (rightEl && inchi2) _renderCard(rightEl, inchi2);
}

function _renderCard(el, inchi) {
    el.innerHTML = "<div class='mol-loading'></div>";

    fetch("http://127.0.0.1:5000/api/render_3d_image", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
            inchi:  inchi,
            width:  Math.max(el.offsetWidth  || 300, 150),
            height: Math.max(el.offsetHeight || 200, 120)
        })
    })
    .then(function(r) { return r.json(); })
    .then(function(data) {
        if (data.error) throw new Error(data.error);
        if (data.sdf) _sdfCache[inchi] = data.sdf;

        el.innerHTML = "";
        var img = document.createElement("img");
        img.src   = "data:image/png;base64," + data.image;
        img.style.cssText = "width:100%;height:100%;object-fit:contain;display:block;cursor:pointer;";
        img.title = "Click for interactive 3D";

        var hint = document.createElement("div");
        hint.className = "mol-hint";
        hint.textContent = "click for 3D \u2197";

        el.appendChild(img);
        el.appendChild(hint);

        img.addEventListener("click", function(e) {
            e.stopPropagation();
            _openModal(inchi);
        });
    })
    .catch(function(err) {
        console.info("PNG failed, 2D fallback:", err.message);
        _render2D(el, inchi);
    });
}

function _render2D(el, inchi) {
    _fetchSDF(inchi).then(function(sdf) {
        if (!sdf) { el.innerHTML = "<div class='mol-error'>Could not render</div>"; return; }
        _rdkitReady.then(function() {
            try {
                var mol = RDKit.get_mol(sdf);
                if (!mol || !mol.is_valid()) throw new Error();
                var svg = mol.get_svg();
                mol.delete();
                el.innerHTML = svg;
                var s = el.querySelector("svg");
                if (s) {
                    s.removeAttribute("width");
                    s.removeAttribute("height");
                    s.style.cssText = "width:100%;height:100%;cursor:pointer;display:block;";
                }
                el.style.cursor = "pointer";
                el.title = "Click for interactive 3D";
                el.addEventListener("click", function(e) { e.stopPropagation(); _openModal(inchi); });
            } catch(e) {
                el.innerHTML = "<div class='mol-error'>Could not render</div>";
            }
        });
    });
}

// ── Modal — built once, reused forever ────────────────
// The 3Dmol viewer is created once inside _modalCanvas and
// NEVER destroyed. Switching molecules only calls clear() +
// addModel() on the existing viewer. This is why multi-open works.

function _buildModal() {
    if (_modalBuilt) return;
    _modalBuilt = true;

    var backdrop = document.createElement("div");
    backdrop.className = "mol-modal-backdrop";
    backdrop.style.display = "none";

    var box = document.createElement("div");
    box.className = "mol-modal";

    var closeBtn = document.createElement("button");
    closeBtn.className = "mol-modal-close";
    closeBtn.innerHTML = "&#x2715;";

    var inchiLabel = document.createElement("div");
    inchiLabel.className = "mol-modal-header";

    var canvasWrap = document.createElement("div");
    canvasWrap.className = "mol-modal-body";

    box.appendChild(closeBtn);
    box.appendChild(inchiLabel);
    box.appendChild(canvasWrap);
    backdrop.appendChild(box);
    document.body.appendChild(backdrop);

    closeBtn.addEventListener("click", function() { backdrop.style.display = "none"; });
    backdrop.addEventListener("click", function(ev) { if (ev.target === backdrop) backdrop.style.display = "none"; });
    document.addEventListener("keydown", function(ev) {
        if (ev.key === "Escape" && backdrop.style.display !== "none") backdrop.style.display = "none";
    });

    _modalBackdrop   = backdrop;
    _modalInchiLabel = inchiLabel;
    _modalCanvas     = canvasWrap;

    // Create the 3Dmol viewer once here — it lives permanently inside canvasWrap
    _modalViewer = $3Dmol.createViewer(canvasWrap, { backgroundColor: "white" });
}

function _openModal(inchi) {
    _buildModal();

    _modalInchiLabel.textContent = inchi;
    _modalBackdrop.style.display = "flex";

    var loadInto = function(sdf) {
        _modalViewer.clear();
        _modalViewer.removeAllModels();
        _modalViewer.addModel(sdf, "sdf");
        _modalViewer.setStyle({}, { stick: {}, sphere: { scale: 0.3 } });
        _modalViewer.zoomTo();
        _modalViewer.render();
        setTimeout(function() {
            try { _modalViewer.resize(); _modalViewer.render(); } catch(e) {}
        }, 80);
    };

    if (_sdfCache[inchi]) {
        loadInto(_sdfCache[inchi]);
        return;
    }

    _modalCanvas.innerHTML = "<div class='mol-loading' style='width:100%;height:100%;display:flex;align-items:center;justify-content:center;'></div>";

    _fetchSDF(inchi).then(function(sdf) {
        if (!sdf) {
            _modalCanvas.innerHTML = "<div class='mol-error' style='padding:20px;text-align:center;'>Could not load 3D structure</div>";
            return;
        }
        _modalCanvas.innerHTML = "";
        _modalViewer = $3Dmol.createViewer(_modalCanvas, { backgroundColor: "white" });
        loadInto(sdf);
    });
}

function _fetchSDF(inchi) {
    if (!inchi) return Promise.resolve(null);
    if (_sdfCache[inchi]) return Promise.resolve(_sdfCache[inchi]);
    if (_fetching[inchi]) return _fetching[inchi];

    var p = fetch("http://127.0.0.1:5000/api/generate_3d", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ inchi: inchi })
    })
    .then(function(r) { return r.json(); })
    .then(function(data) {
        delete _fetching[inchi];
        if (data.error) throw new Error(data.error);
        _sdfCache[inchi] = data.sdf;
        return data.sdf;
    })
    .catch(function() { delete _fetching[inchi]; return null; });

    _fetching[inchi] = p;
    return p;
}