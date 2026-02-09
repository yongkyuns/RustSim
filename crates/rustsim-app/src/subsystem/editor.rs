//! Subsystem editing operations.

use std::collections::HashSet;

use rustsim_types::{
    NodeInstance, PortDirection, Position, SimulationGraph, SubsystemGraph,
};

/// Editor for creating and manipulating subsystems
pub struct SubsystemEditor;

impl SubsystemEditor {
    /// Create a new subsystem node at the given position
    pub fn create_subsystem(
        graph: &mut SimulationGraph,
        id: String,
        name: &str,
        position: Position,
    ) -> String {
        let mut node = NodeInstance::new(id.clone(), "Subsystem".to_string(), position);
        node.name = name.to_string();

        // Initialize with an empty subsystem graph
        node.graph = Some(Box::new(SubsystemGraph::default()));

        graph.add_node(node);
        id
    }

    /// Delete a subsystem and all its contents
    pub fn delete_subsystem(graph: &mut SimulationGraph, subsystem_id: &str) {
        graph.remove_node(subsystem_id);
    }

    /// Add an interface node to a subsystem's internal graph
    ///
    /// Interface nodes define the boundary between the subsystem and its parent.
    /// There should typically be only one interface node per subsystem.
    pub fn add_interface_node(
        subsystem_graph: &mut SubsystemGraph,
        id: String,
        position: Position,
    ) -> String {
        let node = NodeInstance::new(id.clone(), "Interface".to_string(), position);

        subsystem_graph.nodes.push(node);
        id
    }

    /// Add a port to an interface node
    ///
    /// For interface nodes:
    /// - Adding an INPUT creates an output on the parent subsystem (data flows OUT of subsystem)
    /// - Adding an OUTPUT creates an input on the parent subsystem (data flows INTO subsystem)
    pub fn add_port_to_interface(
        interface_node: &mut NodeInstance,
        name: &str,
        direction: PortDirection,
        parent_subsystem: Option<&mut NodeInstance>,
    ) {
        match direction {
            PortDirection::Input => {
                // Input on interface = Output on parent
                interface_node.add_input(name);
                if let Some(parent) = parent_subsystem {
                    parent.add_output(name);
                }
            }
            PortDirection::Output => {
                // Output on interface = Input on parent
                interface_node.add_output(name);
                if let Some(parent) = parent_subsystem {
                    parent.add_input(name);
                }
            }
        }
    }

    /// Remove a port from an interface node
    pub fn remove_port_from_interface(
        interface_node: &mut NodeInstance,
        direction: PortDirection,
        parent_subsystem: Option<&mut NodeInstance>,
    ) -> bool {
        match direction {
            PortDirection::Input => {
                let removed = interface_node.remove_input();
                if removed {
                    if let Some(parent) = parent_subsystem {
                        parent.remove_output();
                    }
                }
                removed
            }
            PortDirection::Output => {
                let removed = interface_node.remove_output();
                if removed {
                    if let Some(parent) = parent_subsystem {
                        parent.remove_input();
                    }
                }
                removed
            }
        }
    }

    /// Convert a selection of nodes into a subsystem
    ///
    /// This creates a new subsystem node and moves the selected nodes into it.
    /// Connections between selected nodes are moved into the subsystem.
    /// Connections crossing the boundary become interface connections.
    pub fn convert_selection_to_subsystem(
        graph: &mut SimulationGraph,
        selected_ids: &[String],
        subsystem_id: String,
        subsystem_name: &str,
    ) -> Result<String, String> {
        if selected_ids.is_empty() {
            return Err("No nodes selected".to_string());
        }

        // Calculate the center position of selected nodes
        let mut sum_x = 0.0;
        let mut sum_y = 0.0;
        let mut count = 0;

        for id in selected_ids {
            if let Some(node) = graph.get_node(id) {
                sum_x += node.position.x;
                sum_y += node.position.y;
                count += 1;
            }
        }

        let center = if count > 0 {
            Position::new(sum_x / count as f32, sum_y / count as f32)
        } else {
            Position::zero()
        };

        // Create the subsystem node
        let mut subsystem_node = NodeInstance::new(
            subsystem_id.clone(),
            "Subsystem".to_string(),
            center,
        );
        subsystem_node.name = subsystem_name.to_string();

        // Create internal graph
        let mut internal_graph = SubsystemGraph::default();
        let selected_set: HashSet<String> = selected_ids.iter().cloned().collect();

        // Move selected nodes into the subsystem
        for id in selected_ids {
            if let Some(node) = graph.nodes.remove(id) {
                // Adjust position relative to center
                let mut adjusted_node = node;
                adjusted_node.position.x -= center.x;
                adjusted_node.position.y -= center.y;

                internal_graph.nodes.push(adjusted_node);
            }
        }

        // Move internal connections
        let mut external_connections = Vec::new();
        graph.connections.retain(|conn| {
            let source_internal = selected_set.contains(&conn.source_node_id);
            let target_internal = selected_set.contains(&conn.target_node_id);

            if source_internal && target_internal {
                // Both nodes are inside - move to internal graph
                internal_graph.connections.push(conn.clone());
                false // Remove from parent graph
            } else if source_internal || target_internal {
                // Boundary connection - needs interface handling
                external_connections.push(conn.clone());
                false // Remove from parent graph for now
            } else {
                // External connection - keep in parent
                true
            }
        });

        // TODO: Handle boundary connections by creating interface node and ports
        // For now, we just lose the external connections (simplified version)

        subsystem_node.graph = Some(Box::new(internal_graph));
        graph.add_node(subsystem_node);

        Ok(subsystem_id)
    }

    /// Expand a subsystem back into individual nodes
    ///
    /// This is the inverse of convert_selection_to_subsystem.
    pub fn expand_subsystem(
        graph: &mut SimulationGraph,
        subsystem_id: &str,
    ) -> Result<Vec<String>, String> {
        let subsystem = graph
            .nodes
            .remove(subsystem_id)
            .ok_or_else(|| format!("Subsystem '{}' not found", subsystem_id))?;

        if subsystem.block_type != "Subsystem" {
            return Err(format!("Node '{}' is not a subsystem", subsystem_id));
        }

        let internal_graph = subsystem
            .graph
            .ok_or_else(|| "Subsystem has no internal graph".to_string())?;

        let center = subsystem.position;
        let mut expanded_ids = Vec::new();

        // Move nodes from subsystem back to parent
        for mut node in internal_graph.nodes {
            // Adjust position back to parent coordinates
            node.position.x += center.x;
            node.position.y += center.y;

            expanded_ids.push(node.id.clone());
            graph.add_node(node);
        }

        // Move connections back
        for conn in internal_graph.connections {
            graph.add_connection(conn);
        }

        // Remove connections to/from the subsystem node
        graph.connections.retain(|conn| {
            conn.source_node_id != subsystem_id && conn.target_node_id != subsystem_id
        });

        Ok(expanded_ids)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_subsystem() {
        let mut graph = SimulationGraph::new();

        let id = SubsystemEditor::create_subsystem(
            &mut graph,
            "sub1".to_string(),
            "Test Subsystem",
            Position::new(100.0, 100.0),
        );

        assert_eq!(id, "sub1");
        assert!(graph.get_node("sub1").is_some());

        let node = graph.get_node("sub1").unwrap();
        assert_eq!(node.block_type, "Subsystem");
        assert_eq!(node.name, "Test Subsystem");
        assert!(node.graph.is_some());
    }

    #[test]
    fn test_delete_subsystem() {
        let mut graph = SimulationGraph::new();
        SubsystemEditor::create_subsystem(
            &mut graph,
            "sub1".to_string(),
            "Test",
            Position::zero(),
        );

        assert!(graph.get_node("sub1").is_some());

        SubsystemEditor::delete_subsystem(&mut graph, "sub1");
        assert!(graph.get_node("sub1").is_none());
    }

    #[test]
    fn test_add_interface_node() {
        let mut subsystem_graph = SubsystemGraph::default();

        let id = SubsystemEditor::add_interface_node(
            &mut subsystem_graph,
            "if1".to_string(),
            Position::zero(),
        );

        assert_eq!(id, "if1");
        assert_eq!(subsystem_graph.nodes.len(), 1);
        assert_eq!(subsystem_graph.nodes[0].block_type, "Interface");
    }

    #[test]
    fn test_add_port_to_interface() {
        let mut interface = NodeInstance::new(
            "if1".to_string(),
            "Interface".to_string(),
            Position::zero(),
        );
        let mut parent = NodeInstance::new(
            "sub1".to_string(),
            "Subsystem".to_string(),
            Position::zero(),
        );

        // Add input to interface (output on parent)
        SubsystemEditor::add_port_to_interface(
            &mut interface,
            "in1",
            PortDirection::Input,
            Some(&mut parent),
        );

        assert_eq!(interface.inputs.len(), 1);
        assert_eq!(parent.outputs.len(), 1);

        // Add output to interface (input on parent)
        SubsystemEditor::add_port_to_interface(
            &mut interface,
            "out1",
            PortDirection::Output,
            Some(&mut parent),
        );

        assert_eq!(interface.outputs.len(), 1);
        assert_eq!(parent.inputs.len(), 1);
    }

    #[test]
    fn test_convert_selection_to_subsystem() {
        let mut graph = SimulationGraph::new();

        // Add some nodes
        let node1 = NodeInstance::new(
            "node1".to_string(),
            "Constant".to_string(),
            Position::new(0.0, 0.0),
        );
        let node2 = NodeInstance::new(
            "node2".to_string(),
            "Constant".to_string(),
            Position::new(100.0, 100.0),
        );

        graph.add_node(node1);
        graph.add_node(node2);

        let selected = vec!["node1".to_string(), "node2".to_string()];

        let result = SubsystemEditor::convert_selection_to_subsystem(
            &mut graph,
            &selected,
            "subsys1".to_string(),
            "New Subsystem",
        );

        assert!(result.is_ok());

        // Original nodes should be removed
        assert!(graph.get_node("node1").is_none());
        assert!(graph.get_node("node2").is_none());

        // Subsystem should exist
        let subsys = graph.get_node("subsys1").unwrap();
        assert_eq!(subsys.block_type, "Subsystem");

        // Nodes should be in subsystem
        let internal = subsys.graph.as_ref().unwrap();
        assert_eq!(internal.nodes.len(), 2);
    }

    #[test]
    fn test_expand_subsystem() {
        let mut graph = SimulationGraph::new();

        // Create subsystem with nodes
        let mut subsystem = NodeInstance::new(
            "sub1".to_string(),
            "Subsystem".to_string(),
            Position::new(100.0, 100.0),
        );

        let mut internal = SubsystemGraph::default();
        internal.nodes.push(NodeInstance::new(
            "node1".to_string(),
            "Constant".to_string(),
            Position::new(0.0, 0.0),
        ));
        internal.nodes.push(NodeInstance::new(
            "node2".to_string(),
            "Constant".to_string(),
            Position::new(50.0, 50.0),
        ));

        subsystem.graph = Some(Box::new(internal));
        graph.add_node(subsystem);

        let result = SubsystemEditor::expand_subsystem(&mut graph, "sub1");
        assert!(result.is_ok());

        let expanded = result.unwrap();
        assert_eq!(expanded.len(), 2);

        // Subsystem should be removed
        assert!(graph.get_node("sub1").is_none());

        // Nodes should be back in parent
        assert!(graph.get_node("node1").is_some());
        assert!(graph.get_node("node2").is_some());

        // Check positions were adjusted
        let node1 = graph.get_node("node1").unwrap();
        assert_eq!(node1.position.x, 100.0); // 0 + 100
        assert_eq!(node1.position.y, 100.0); // 0 + 100
    }
}
