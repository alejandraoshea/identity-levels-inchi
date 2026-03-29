from flask import Flask, request, jsonify
from backend.inchi.compare import compare_files
from backend.inchi.determine_levels_id import InChi
from backend.inchi.config_loader import load_config
from backend.routes.inchi_comparison_routes import inchi_comparison_routes
import tempfile
import json
from flask_cors import CORS

app = Flask(__name__)
CORS(app)
CORS(app, resources={r"/*": {"origins": "http://localhost:3000"}})

app.register_blueprint(inchi_comparison_routes)

if __name__ == "__main__":
    app.run(host="127.0.0.1", port=8080)