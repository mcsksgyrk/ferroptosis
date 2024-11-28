// Initialize the graph visualization
function initGraph() {
    const container = document.getElementById('container');
    const width = container.clientWidth;
    const height = container.clientHeight;

    const svg = d3.select("svg")
        .attr("width", "100%")
        .attr("height", "100%")
        .attr("viewBox", [0, 0, width, height]);

    const zoom = d3.zoom()
        .scaleExtent([0.1, 4])
        .on("zoom", (event) => {
            g.attr("transform", event.transform);
        });

    svg.call(zoom);
    const g = svg.append("g");

    return { g, width, height };
}

// Create force simulation
function createSimulation(nodes, edges, width, height) {
    return d3.forceSimulation(nodes)
        .force("link", d3.forceLink(edges)
            .id(d => d.id)
            .distance(100))
        .force("charge", d3.forceManyBody()
            .strength(-1000))
        .force("center", d3.forceCenter(width / 2, height / 2))
        .force("x", d3.forceX(width / 2).strength(0.1))
        .force("y", d3.forceY(height / 2).strength(0.1))
        .force("collision", d3.forceCollide().radius(50));
}

// Highlight connected nodes and links
function highlightConnections(event, d, node, link, data) {
    // Fade all nodes and links first
    node.classed("faded", true);
    link.classed("faded", true);

    // Find connected elements
    const connectedLinks = data.edges.filter(l =>
        l.source.id === d.id || l.target.id === d.id
    );
    const connectedNodes = new Set();
    connectedLinks.forEach(l => {
        connectedNodes.add(l.source.id);
        connectedNodes.add(l.target.id);
    });

    // Highlight connected nodes
    node.filter(n => connectedNodes.has(n.id))
        .classed("faded", false)
        .select("circle")
        .attr("fill", "#FF6B6B")
        .attr("r", 25)
        .attr("stroke", "#FF4949")
        .attr("stroke-width", 4);

    // Highlight the selected node specially
    node.filter(n => n.id === d.id)
        .select("circle")
        .attr("fill", "#FFD93D")
        .attr("r", 30)
        .attr("stroke", "#FF8C00")
        .attr("stroke-width", 5);

    // Highlight connected node labels
    node.filter(n => connectedNodes.has(n.id))
        .select("text")
        .attr("font-size", "16px")
        .attr("font-weight", "bold")
        .attr("fill", "#FF4949");

    // Highlight connected links
    connectedLinks.forEach(l => {
        d3.select(`#link-${l.source.id}-${l.target.id}`)
            .classed("faded", false)
            .attr("stroke", "#FF4949")
            .attr("stroke-width", 5)
            .attr("stroke-opacity", 1);
    });
}

// Reset highlighting
function resetHighlight(node, link) {
    node.classed("faded", false)
        .select("circle")
        .attr("fill", "#69b3a2")
        .attr("r", 20)
        .attr("stroke", "#333")
        .attr("stroke-width", 3);

    node.select("text")
        .attr("font-size", "14px")
        .attr("font-weight", "normal")
        .attr("fill", "black");

    link.classed("faded", false)
        .attr("stroke", "#999")
        .attr("stroke-width", 3)
        .attr("stroke-opacity", 0.6);
}

// Drag functions
function dragstarted(event, d, simulation) {
    if (!event.active) simulation.alphaTarget(0.3).restart();
    d.fx = d.x;
    d.fy = d.y;
}

function dragged(event, d) {
    d.fx = event.x;
    d.fy = event.y;
}

function dragended(event, d, simulation) {
    if (!event.active) simulation.alphaTarget(0);
    d.fx = null;
    d.fy = null;
}

// Main function to create the graph
function createGraph(data) {
    const { g, width, height } = initGraph();
    const simulation = createSimulation(data.nodes, data.edges, width, height);

    // Create links
    const link = g.selectAll(".link")
        .data(data.edges)
        .join("line")
        .attr("class", "link")
        .attr("id", d => `link-${d.source.id}-${d.target.id}`)
        .attr("stroke-width", 3);

    // Create nodes
    const node = g.selectAll(".node")
        .data(data.nodes)
        .join("g")
        .attr("class", "node")
        .call(d3.drag()
            .on("start", (event, d) => dragstarted(event, d, simulation))
            .on("drag", dragged)
            .on("end", (event, d) => dragended(event, d, simulation)));

    // Add circles to nodes
    node.append("circle")
        .attr("r", 20)
        .attr("stroke-width", 3);

    // Add labels to nodes
    node.append("text")
        .text(d => d.name)
        .attr("x", 25)
        .attr("y", 5)
        .attr("font-size", "14px");

    // Add event listeners
    node.on("mouseover", (event, d) => highlightConnections(event, d, node, link, data))
        .on("mouseout", () => resetHighlight(node, link));

    // Update positions on each tick
    simulation.on("tick", () => {
        link
            .attr("x1", d => d.source.x)
            .attr("y1", d => d.source.y)
            .attr("x2", d => d.target.x)
            .attr("y2", d => d.target.y);

        node.attr("transform", d => `translate(${d.x},${d.y})`);
    });
}

// Fetch data and initialize the graph
fetch('/graph-data')
    .then(response => response.json())
    .then(data => createGraph(data));
