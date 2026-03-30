from flask import Blueprint, jsonify, request
from backend.inchi.determine_levels_id import InChi
from backend.inchi.compare import compare_files
from backend.inchi.config_loader import load_config, build_config_from_levels
import tempfile
from dotenv import load_dotenv

load_dotenv()

inchi_comparison_routes = Blueprint("inchi_comparison_routes", __name__)

@inchi_comparison_routes.route("/api/compare_inchis", methods=["POST"])
def compare_inchis():
    try:
        data = request.get_json()

        if not data:
            return jsonify({"message": "No JSON received"}), 400

        inchi1 = data.get("inchi1")
        inchi2 = data.get("inchi2")

        if not inchi1 or not inchi2:
            return jsonify({"message": "Missing InChIs"}), 400

        config = load_config()

        comparison = InChi.get_ids(inchi1, inchi2, config)

        results = {k.name: v for k, v in comparison.items()}

        return jsonify({"results": results})

    except Exception as e:
        print("ERROR:", e)
        return jsonify({"message": str(e)}), 500
  
@inchi_comparison_routes.route("/api/compare_inchis_custom", methods=["POST"])
def compare_inchis_custom():
    try:
        data = request.get_json()

        if not data:
            return jsonify({"message": "No JSON received"}), 400

        inchi1 = data.get("inchi1")
        inchi2 = data.get("inchi2")
        selected_levels = data.get("levels", [])

        if not inchi1 or not inchi2:
            return jsonify({"message": "Missing InChIs"}), 400

        base_config = load_config()
        config = build_config_from_levels(selected_levels, base_config)

        comparison = InChi.get_ids(inchi1, inchi2, config)

        results = {k.name: v for k, v in comparison.items()}

        return jsonify({"results": results})

    except Exception as e:
        print("ERROR:", e)
        return jsonify({"message": str(e)}), 500


@inchi_comparison_routes.route("/api/inchi_levels", methods=["GET"])
def get_inchi_levels():
    return jsonify([
        {"key": "complete_identity", "label": "Complete Identity"},
        {"key": "isotope", "label": "Isotope Independence"},
        {"key": "salt", "label": "Salt Independence"},
        {"key": "charge", "label": "Charge Independence"},
        {"key": "double_bond", "label": "Double Bond Independence"},
        {"key": "cis_trans", "label": "Cis/Trans Independence"},
        {"key": "tautomer", "label": "Tautomer Independence"},
        {"key": "substituent", "label": "Substituent Independence"}
    ])

@inchi_comparison_routes.route("/api/compare_files", methods=["POST"])
def compare_files_api():
    try:
        data = request.get_json()

        list1 = data.get("list1", [])
        list2 = data.get("list2", [])

        if not list1 or not list2:
            return jsonify({"message": "Both lists required"}), 400

        config = load_config()

        results = []

        for i1 in list1:
            for i2 in list2:
                comparison = InChi.get_ids(i1, i2, config)

                results.append({
                    "inchi_1": i1,
                    "inchi_2": i2,
                    "results": {k.name: v for k, v in comparison.items()}
                })

        return jsonify({"comparisons": results})

    except Exception as e:
        print("ERROR:", e)
        return jsonify({"message": str(e)}), 500