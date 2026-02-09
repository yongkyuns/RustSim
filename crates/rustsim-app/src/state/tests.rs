//! Tests for state modules

#[cfg(test)]
mod tests {
    use super::super::*;
    use rustsim_types::{Position, SimulationGraph, SimulationSettings};

    #[test]
    fn test_graph_state_creation() {
        let graph_state = GraphState::new();
        assert!(graph_state.graph.nodes.is_empty());
        assert!(graph_state.graph.connections.is_empty());
        assert!(!graph_state.block_types.is_empty());
    }

    #[test]
    fn test_graph_state_add_node() {
        let mut graph_state = GraphState::new();
        let node_id = graph_state.add_node("Constant", Position::new(100.0, 100.0));

        assert!(!node_id.is_empty());
        assert_eq!(graph_state.graph.nodes.len(), 1);
        assert!(graph_state.graph.get_node(&node_id).is_some());
    }

    #[test]
    fn test_graph_state_add_connection() {
        let mut graph_state = GraphState::new();
        let node1 = graph_state.add_node("Constant", Position::new(100.0, 100.0));
        let node2 = graph_state.add_node("Amplifier", Position::new(200.0, 100.0));

        let conn_id = graph_state.add_connection(&node1, 0, &node2, 0);

        assert!(!conn_id.is_empty());
        assert_eq!(graph_state.graph.connections.len(), 1);
    }

    #[test]
    fn test_graph_state_port_management() {
        let mut graph_state = GraphState::new();
        let node_id = graph_state.add_node("Adder", Position::new(100.0, 100.0));

        // Adder starts with 2 inputs
        let node = graph_state.graph.get_node(&node_id).unwrap();
        assert_eq!(node.inputs.len(), 2);

        // Should be able to add more inputs (Adder has unlimited inputs)
        assert!(graph_state.can_add_input(&node_id));
        assert!(graph_state.add_input_port(&node_id));

        let node = graph_state.graph.get_node(&node_id).unwrap();
        assert_eq!(node.inputs.len(), 3);
    }

    #[test]
    fn test_simulation_state_lifecycle() {
        let mut sim_state = SimulationState::new();

        assert!(!sim_state.is_running());
        assert_eq!(sim_state.sim_time, 0.0);

        sim_state.start();
        assert!(sim_state.is_running());

        sim_state.stop();
        assert!(!sim_state.is_running());

        sim_state.reset();
        assert_eq!(sim_state.sim_time, 0.0);
    }

    #[test]
    fn test_simulation_state_execution_order() {
        let mut graph = SimulationGraph::new();
        let mut sim_state = SimulationState::new();

        // Create a simple chain: Constant -> Amplifier -> Scope
        graph.add_node(rustsim_types::NodeInstance::new(
            "const1".to_string(),
            "Constant".to_string(),
            Position::new(0.0, 0.0),
        ));

        graph.add_node(rustsim_types::NodeInstance::new(
            "amp1".to_string(),
            "Amplifier".to_string(),
            Position::new(100.0, 0.0),
        ));

        graph.add_node(rustsim_types::NodeInstance::new(
            "scope1".to_string(),
            "Scope".to_string(),
            Position::new(200.0, 0.0),
        ));

        // Add connections
        graph.add_connection(rustsim_types::Connection::new(
            "conn1".to_string(),
            "const1".to_string(),
            0,
            "amp1".to_string(),
            0,
        ));

        graph.add_connection(rustsim_types::Connection::new(
            "conn2".to_string(),
            "amp1".to_string(),
            0,
            "scope1".to_string(),
            0,
        ));

        sim_state.create_blocks_from_graph(&graph);
        let order = sim_state.execution_order();

        // Ensure topological order: const1 before amp1 before scope1
        let const_pos = order.iter().position(|id| id == "const1").unwrap();
        let amp_pos = order.iter().position(|id| id == "amp1").unwrap();
        let scope_pos = order.iter().position(|id| id == "scope1").unwrap();

        assert!(const_pos < amp_pos);
        assert!(amp_pos < scope_pos);
    }

    #[test]
    fn test_canvas_state_selection() {
        let mut canvas_state = CanvasState::new();

        assert!(canvas_state.selected_nodes.is_empty());

        canvas_state.select_node("node1".to_string());
        assert_eq!(canvas_state.selected_nodes.len(), 1);
        assert!(canvas_state.is_node_selected("node1"));

        canvas_state.deselect_node("node1");
        assert!(canvas_state.selected_nodes.is_empty());
    }

    #[test]
    fn test_canvas_state_zoom() {
        let mut canvas_state = CanvasState::new();
        let initial_zoom = canvas_state.zoom;

        canvas_state.zoom_in();
        assert!(canvas_state.zoom > initial_zoom);

        canvas_state.zoom_out();
        assert!((canvas_state.zoom - initial_zoom).abs() < 0.01);
    }

    #[test]
    fn test_plot_state_data_collection() {
        let mut plot_state = PlotState::new();

        assert!(plot_state.plot_data.is_empty());

        plot_state.add_plot_point(0.0, vec![1.0, 2.0, 3.0]);
        plot_state.add_plot_point(0.1, vec![1.1, 2.1, 3.1]);

        assert_eq!(plot_state.plot_data.len(), 2);
        assert_eq!(plot_state.plot_data[0].0, 0.0);
        assert_eq!(plot_state.plot_data[0].1, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_plot_state_history() {
        let mut plot_state = PlotState::new();

        plot_state.add_plot_point(0.0, vec![1.0]);
        plot_state.add_plot_point(0.1, vec![1.1]);

        plot_state.save_to_history();
        assert_eq!(plot_state.result_history.len(), 1);

        // Clear and add new data
        plot_state.clear();
        plot_state.add_plot_point(0.0, vec![2.0]);
        plot_state.save_to_history();

        assert_eq!(plot_state.result_history.len(), 2);
    }

    #[test]
    fn test_history_state_undo_redo() {
        let mut history_state = HistoryState::new();
        let graph1 = SimulationGraph::new();
        let settings1 = SimulationSettings::default();

        // Save initial state
        history_state.save_state(&graph1, &settings1);

        // Create modified state
        let mut graph2 = SimulationGraph::new();
        graph2.add_node(rustsim_types::NodeInstance::new(
            "node1".to_string(),
            "Constant".to_string(),
            Position::new(0.0, 0.0),
        ));

        // Undo should restore graph1
        let (restored_graph, _) = history_state.undo(&graph2, &settings1).unwrap();
        assert_eq!(restored_graph.nodes.len(), 0);

        // Redo should restore graph2
        let (restored_graph, _) = history_state.redo(&graph1, &settings1).unwrap();
        assert_eq!(restored_graph.nodes.len(), 1);
    }

    #[test]
    fn test_app_state_composition() {
        let app_state = AppState::new();

        // All sub-states should be initialized
        assert!(!app_state.graph_state.block_types.is_empty());
        assert_eq!(app_state.canvas_state.zoom, 1.0);
        assert!(!app_state.simulation_state.is_running());
        assert!(app_state.plot_state.plot_data.is_empty());
        assert_eq!(app_state.history_state.undo_size(), 0);
    }

    #[test]
    fn test_app_state_delegate_methods() {
        let app_state = AppState::new();

        // Test graph access
        let graph = app_state.graph();
        assert!(graph.nodes.is_empty());

        // Test block type access
        let const_type = app_state.get_block_type("Constant");
        assert!(const_type.is_some());
        assert_eq!(const_type.unwrap().name, "Constant");
    }
}
