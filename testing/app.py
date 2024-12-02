# app.py
from flask import Flask, jsonify, request
from flask_cors import CORS
import sqlite3

app = Flask(__name__)
CORS(app)

def get_db_connection():
    conn = sqlite3.connect('../../output/omnipath.db')
    conn.row_factory = sqlite3.Row
    return conn

@app.route('/api/nodes')
def get_nodes():
    conn = get_db_connection()
    cursor = conn.cursor()

    # Support filtering by parameters
    filters = request.args.get('filters', '')
    if filters:
        query = "SELECT * FROM nodes WHERE " + filters
    else:
        query = "SELECT * FROM nodes"

    cursor.execute(query)
    nodes = [dict(row) for row in cursor.fetchall()]

    conn.close()
    return jsonify(nodes)

@app.route('/api/edges')
def get_edges():
    conn = get_db_connection()
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM edges")
    edges = [dict(row) for row in cursor.fetchall()]
    conn.close()
    return jsonify(edges)

@app.route('/api/node/<node_id>')
def get_node_info(node_id):
    conn = get_db_connection()
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM nodes WHERE id = ?", (node_id,))
    node = dict(cursor.fetchone())
    conn.close()
    return jsonify(node)

if __name__ == '__main__':
    app.run(debug=True)
