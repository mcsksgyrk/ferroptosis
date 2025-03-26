DROP TABLE IF EXISTS `node`;
DROP TABLE IF EXISTS `edge`;
DROP TABLE IF EXISTS `node_identifier`;
DROP TABLE IF EXISTS `cellline`;
DROP TABLE IF EXISTS `edge_cellline_map`;
DROP TABLE IF EXISTS `disease`;
DROP TABLE IF EXISTS `disease_edge`;
PRAGMA foreign_keys = ON;

CREATE TABLE `node` (
    `id` INTEGER PRIMARY KEY,
    `name` TEXT NOT NULL,         -- Primary identifier (UniProt/PubChem/KEGG ID)
    `primary_id_type` TEXT,       -- Type of the primary identifier in name field
    `display_name` TEXT,          -- Human-readable name
    `tax_id` INTEGER,             -- Only for biological entities, now nullable
    `type` TEXT NOT NULL DEFAULT 'protein', -- 'protein', 'small_molecule', 'rna', etc.
    `pathways` TEXT,
    `source` TEXT,
    `function` TEXT
);

CREATE TABLE `node_identifier` (
    `node_id` INTEGER NOT NULL,
    `id_type` TEXT NOT NULL,      -- 'kegg_id', 'uniprot_id', 'pubchem_id', etc.
    `is_primary` BOOLEAN DEFAULT 0, -- Indicates if this is the primary ID used in node.name
    `id_value` TEXT NOT NULL,
    PRIMARY KEY (`node_id`, `id_type`),
    FOREIGN KEY (`node_id`) REFERENCES `node`(`id`) ON UPDATE NO ACTION ON DELETE CASCADE
);

CREATE TABLE `cellline` (
    `id` INTEGER PRIMARY KEY,
    `cellline_id` TEXT NOT NULL UNIQUE,  -- Original identifier from source database
    `cellline_name` TEXT,         -- Human-readable name
    `cellline_species` TEXT,      -- Species of origin (e.g., human, mouse)
    `cellline_disease` TEXT,      -- Associated disease if applicable
    `cell_describe` TEXT          -- Additional information about the cell line
);

CREATE TABLE `disease` (
    `id` INTEGER PRIMARY KEY,
    `disease_id` TEXT NOT NULL UNIQUE,   -- External identifier (e.g., OMIM, DOID)
    `disease_name` TEXT NOT NULL,        -- Human-readable name
    `disease_type` TEXT,                 -- Classification or category
    `description` TEXT
);

CREATE TABLE `edge` (
    `id` INTEGER PRIMARY KEY,
    `interactor_a_node_id` INTEGER NOT NULL,
    `interactor_b_node_id` INTEGER NOT NULL,
    `interactor_a_node_name` TEXT NOT NULL,
    `interactor_b_node_name` TEXT NOT NULL,
    `layer` TEXT NOT NULL,
    `interaction_types` TEXT,
    `source_db` TEXT NOT NULL,
    FOREIGN KEY(`interactor_a_node_id`) REFERENCES node (`id`) ON UPDATE NO ACTION ON DELETE CASCADE,
    FOREIGN KEY(`interactor_b_node_id`) REFERENCES node (`id`) ON UPDATE NO ACTION ON DELETE CASCADE
);

-- Junction table to connect edges with cell lines (many-to-many relationship)
CREATE TABLE `edge_cellline_map` (
    `edge_id` INTEGER NOT NULL,
    `cellline_id` INTEGER NOT NULL,
    `detection_method` TEXT,       -- Method used to detect the interaction in this cell line
    `source_reference` TEXT,       -- Reference to the source of this cell line annotation
    PRIMARY KEY (`edge_id`, `cellline_id`),
    FOREIGN KEY (`edge_id`) REFERENCES `edge`(`id`) ON UPDATE NO ACTION ON DELETE CASCADE,
    FOREIGN KEY (`cellline_id`) REFERENCES `cellline`(`id`) ON UPDATE NO ACTION ON DELETE CASCADE
);

-- Table for disease-entity relationships
CREATE TABLE `disease_edge` (
    `id` INTEGER PRIMARY KEY,
    `disease_id` INTEGER NOT NULL,   -- References disease.id
    `node_id` INTEGER NOT NULL,      -- References node.id (protein, compound, etc.)
    `node_name` TEXT NOT NULL,       -- Name of the connected entity
    `association_type` TEXT,         -- Nature of the association (causal, risk factor, etc.)
    `score` REAL,                    -- Statistical measure of association strength
    `evidence_type` TEXT,            -- Type of evidence (GWAS, clinical, etc.)
    `source_db` TEXT NOT NULL,       -- Data source
    FOREIGN KEY(`disease_id`) REFERENCES disease (`id`) ON UPDATE NO ACTION ON DELETE CASCADE,
    FOREIGN KEY(`node_id`) REFERENCES node (`id`) ON UPDATE NO ACTION ON DELETE CASCADE
);

-- Add indexes to improve query performance
CREATE INDEX `idx_node_name` ON `node`(`name`);
CREATE INDEX `idx_node_name_tax` ON `node`(`name`, `tax_id`);
CREATE INDEX `idx_node_type` ON `node`(`type`);
CREATE INDEX `idx_node_identifier_type` ON `node_identifier`(`id_type`);
CREATE INDEX `idx_node_identifier_value` ON `node_identifier`(`id_value`);
CREATE INDEX `idx_node_identifier_node` ON `node_identifier`(`node_id`);
CREATE INDEX `idx_edge_interactors` ON `edge`(`interactor_a_node_id`, `interactor_b_node_id`);
CREATE INDEX `idx_edge_layer` ON `edge`(`layer`);
CREATE INDEX `idx_cellline_id` ON `cellline`(`cellline_id`);
CREATE INDEX `idx_cellline_name` ON `cellline`(`cellline_name`);
CREATE INDEX `idx_edge_cellline_edge` ON `edge_cellline_map`(`edge_id`);
CREATE INDEX `idx_edge_cellline_cellline` ON `edge_cellline_map`(`cellline_id`);
CREATE INDEX `idx_disease_id` ON `disease`(`disease_id`);
CREATE INDEX `idx_disease_name` ON `disease`(`disease_name`);
CREATE INDEX `idx_disease_edge_disease` ON `disease_edge`(`disease_id`);
CREATE INDEX `idx_disease_edge_node` ON `disease_edge`(`node_id`);
