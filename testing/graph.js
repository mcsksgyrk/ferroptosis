// graph.js
class GraphVisualization {
    constructor(containerId) {
        this.svg = d3.select(containerId);
        this.width = +this.svg.attr("width");
        this.height = +this.svg.attr("height");

        this.simulation = d3.forceSimulation()
            .force("link", d3.forceLink().id(d => d.id))
            .force("charge", d3.forceManyBody())
            .force("center", d3.forceCenter(this.width / 2, this.height / 2));

        this.nodeGroup = this.svg.append("g");
        this.linkGroup = this.svg.append("g");

        this.initializeZoom();
    }

    initializeZoom() {
        const zoom = d3.zoom()
            .scaleExtent([0.1, 4])
            .on("zoom", (event) => {
                this.nodeGroup.attr("transform", event.transform);
                this.linkGroup.attr("transform", event.transform);
            });

        this.svg.call(zoom);
    }

    async loadData(filters = '') {
        try {
            const [nodesResponse, edgesResponse] = await Promise.all([
                fetch(`/api/nodes?filters=${filters}`),
                fetch('/api/edges')
            ]);

            const nodes = await nodesResponse.json();
            const edges = await edgesResponse.json();

            this.updateVisualization(nodes, edges);
        } catch (error) {
            console.error("Error loading data:", error);
        }
    }

    updateVisualization(nodes, edges) {
        // Update links
        const link = this.linkGroup
            .selectAll("line")
            .data(edges)
            .join("line")
            .attr("stroke", "#999")
            .attr("stroke-opacity", 0.6);

        // Update nodes
        const node = this.nodeGroup
            .selectAll("circle")
            .data(nodes)
            .join("circle")
            .attr("r", 5)
            .attr("fill", "#69b3a2")
            .call(this.drag());

        // Add node interaction
        node.on("click", (event, d) => this.showNodeDetails(d));

        // Update simulation
        this.simulation
            .nodes(nodes)
            .force("link").links(edges);

        this.simulation.on("tick", () => {
            link
                .attr("x1", d => d.source.x)
                .attr("y1", d => d.source.y)
                .attr("x2", d => d.target.x)
                .attr("y2", d => d.target.y);

            node
                .attr("cx", d => d.x)
                .attr("cy", d => d.y);
        });
    }

    drag() {
        return d3.drag()
            .on("start", (event, d) => {
                if (!event.active) this.simulation.alphaTarget(0.3).restart();
                d.fx = d.x;
                d.fy = d.y;
            })
            .on("drag", (event, d) => {
                d.fx = event.x;
                d.fy = event.y;
            })
            .on("end", (event, d) => {
                if (!event.active) this.simulation.alphaTarget(0);
                d.fx = null;
                d.fy = null;
            });
    }

    async showNodeDetails(node) {
        try {
            const response = await fetch(`/api/node/${node.id}`);
            const details = await response.json();

            // Update the details panel (implement this based on your UI needs)
            document.getElementById('nodeDetails').innerHTML = `
                <h3>Node Details</h3>
                <pre>${JSON.stringify(details, null, 2)}</pre>
            `;
        } catch (error) {
            console.error("Error fetching node details:", error);
        }
    }

    // Add methods for filtering and searching
    filterNodes(criteria) {
        this.loadData(encodeURIComponent(criteria));
    }
}

// Initialize and use
document.addEventListener('DOMContentLoaded', () => {
    const graph = new GraphVisualization('#graph');
    graph.loadData();

    // Example of adding search functionality
    document.querySelector('#searchForm').addEventListener('submit', (e) => {
        e.preventDefault();
        const searchTerm = document.querySelector('#searchInput').value;
        graph.filterNodes(`name LIKE '%${searchTerm}%'`);
    });
});
