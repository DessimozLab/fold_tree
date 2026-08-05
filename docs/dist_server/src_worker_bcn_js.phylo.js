/*
 * ATTENTION: The "eval" devtool has been used (maybe by default in mode: "development").
 * This devtool is neither made for production nor for readable output files.
 * It uses "eval()" calls to create a separate source file in the browser devtools.
 * If you are trying to read the output file, select a different devtool (https://webpack.js.org/configuration/devtool/)
 * or disable the default devtool with "devtool: false".
 * If you are looking for production-ready output files, see mode: "production" (https://webpack.js.org/configuration/mode/).
 */
var PhyloIO;
/******/ (() => { // webpackBootstrap
/******/ 	var __webpack_modules__ = ({

/***/ "./src/model.js":
/*!**********************!*\
  !*** ./src/model.js ***!
  \**********************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
eval("__webpack_require__.r(__webpack_exports__);\n/* harmony export */ __webpack_require__.d(__webpack_exports__, {\n/* harmony export */   \"default\": () => (/* binding */ Model)\n/* harmony export */ });\n/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! d3 */ \"./node_modules/d3/index.js\");\n/* harmony import */ var minhashjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! minhashjs */ \"./node_modules/minhashjs/index.js\");\n/* harmony import */ var biojs_io_newick__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! biojs-io-newick */ \"./node_modules/biojs-io-newick/src/index.js\");\n\n\n\n\nvar uid_model = 0\nvar uid_untitle_counter = 0\n;\nconst { parse_nhx } = __webpack_require__(/*! ./utils.js */ \"./src/utils.js\")\n\nclass Model {\n\n    constructor(data, settings, from_raw_data = true) {\n\n        this.zoom;\n        this.settings = {\n            'uid': null,\n            'domain_extended_data' : {},\n            'extended_data_type' : {'Topology': 'num'},\n            'labels' : {'leaf' : new Set(), 'node':new Set()},\n            'colorlabels' :{'leaf' : new Set(), 'node':new Set([\"Topology\"])},\n            'display_leaves' : true,\n            'mirror': false,\n            'name': null,\n            'first_time_render': true,\n            'data_type' : 'newick',\n            'use_branch_lenght' : true,\n            'show_tooltips' : false,\n            'subsample_label' : true,\n            'display_internal_label' : false,\n            'display_internal_label_left_top' : false,\n            'display_internal_label_left_bottom' : false,\n            'display_duplication' : false,\n            'has_branch_lenght' : true,\n            'has_duplications' : false,\n            'dessimode': false,\n            'multiple_search':false,\n            'show_histogram' : false,\n            'align_tip' : false,\n            'use_meta_for_leaf' : true,\n            'use_meta_for_node' : false,\n            'has_histogram_data' : false,\n            'similarity': [],\n            'style': {\n                'font_size_internal' : 14,\n                'color_accessor' : {'leaf' : null, 'node': \"Topology\"},\n                'color_extent_min': {'leaf' : {}, 'node': {\"Topology\":0}},\n                'color_extent_max':{'leaf' : {}, 'node': {\"Topology\":1}},\n                'number_domain':{'leaf' : '3', 'node': '5'},\n                'color_domain':{'leaf' : ['#253494', '#2C7FB8', '#41B6C4', '#C7E9B4', '#FFFFCC'], 'node': ['#253494', '#2C7FB8', '#41B6C4', '#C7E9B4', '#FFFFCC']},\n                'color_domain_default':{'leaf' : ['#253494', '#2C7FB8', '#41B6C4', '#C7E9B4', '#FFFFCC'], 'node': ['#253494', '#2C7FB8', '#41B6C4', '#C7E9B4', '#FFFFCC']},\n        },\n            'tree': {\n                'node_vertical_size' : 30,\n                'node_horizontal_size' : 40,\n                'node_radius' : 6, // move to style\n                'line_width' : 3,// move to style\n                'font_size':14, // move to style\n                'max_depth' : 0,\n            },\n            'collapse_level': 0,\n            'stack' : {\n                'type': 'genes',//'events',\n                'showHistogramValues' : false,\n                'showHistogramSummaryValue' : true,\n                'legendTxtSize' : 12,\n                'margin' : 8,\n                'xInitialRightMargin' : 45,\n                'stackHeight' : 120,\n                'stackWidth' : 30,\n                'maxStackHeight': 'max', // ratio -> stack height fixed | max -> largest data = stack height\n                'has_support' : false,\n                'only_support' : false,\n\n            },\n        }\n\n        if (settings) {\n\n            for(var key in settings) {\n\n                if (key == 'labels_array_leaf'){\n                    var value = settings[key];\n                    this.settings['labels']['leaf'] = new Set(value);\n                }\n                if (key == 'labels_array_node'){\n                    var value = settings[key];\n                    this.settings['labels']['node'] = new Set(value);\n                }\n                if (key == 'colorlabels_array_leaf'){\n                    var value = settings[key];\n                    this.settings['colorlabels']['leaf'] = new Set(value);\n                }\n\n                if (key == 'colorlabels_array_node'){\n                    var value = settings[key];\n                    this.settings['colorlabels']['node'] = new Set(value);\n                }\n\n                else{\n                    var value = settings[key];\n                    this.settings[key] = value;\n                }\n\n\n            }\n\n        }\n\n        this.settings.name = this.settings.name ? this.settings.name : \"Untitled \" + uid_untitle_counter++\n\n        this.uid = null\n        this.input_data = data;\n        this.leaves = []\n\n\n\n        if (from_raw_data){\n            this.uid = uid_model++;\n            this.settings.uid = this.uid;\n            this.data = this.factory(this.parse());\n        }\n        else{\n            this.uid = settings.uid;\n            this.settings.uid = this.uid;\n            this.data = data\n            data.leaves = this.get_leaves(data)\n            this.traverse(data, function(n,c){\n                n.leaves = this.get_leaves(n)\n            })\n        }\n\n        this.data.root = true;\n        this.data.elementS = {}\n        this.data.elementBCN = {}\n        this.rooted = this.data.children.length !== 3\n        this.big_tree = (this.leaves.length > 500)\n\n        // check that histogram data is present and compute\n        if(this.settings.show_histogram && this.data.evolutionaryEvents) {\n            this.settings.has_histogram_data  = true;\n            this.largestGenome =  0;\n            this.largestEvents = 0; // todo\n\n            this.traverse(this.data , function(n,c){\n\n                let g = n.nr_hogs ? n.nr_hogs : n.nr_proteins\n                if (g > this.largestGenome ) {this.largestGenome = g;}\n\n                if (n.evolutionaryEvents){\n\n                    var ga = n.evolutionaryEvents.gained ? n.evolutionaryEvents.gained : 0\n                    var l = n.evolutionaryEvents.lost ? n.evolutionaryEvents.lost : 0\n                    var d = n.evolutionaryEvents.duplications ? n.evolutionaryEvents.duplications : 0\n\n                    let e = ga + l + d\n\n\n                    if (e > this.largestEvents ) {this.largestEvents = e;}\n\n                }\n\n                if (this.settings.stack.has_support){\n\n                    let g_support = n.nr_hogs_support ? n.nr_hogs_support : n.nr_proteins_support\n                    if (g_support > this.largestGenome_support ) {this.largestGenome_support = g_support;}\n\n                    if (n.evolutionaryEvents_support){\n\n\n                        var ga_support = n.evolutionaryEvents_support.gained ? n.evolutionaryEvents_support.gained : 0\n                        var l_support = n.evolutionaryEvents_support.lost ? n.evolutionaryEvents_support.lost : 0\n                        var d_support = n.evolutionaryEvents_support.duplications ? n.evolutionaryEvents_support.duplications : 0\n\n                        let e_support = ga_support + l_support + d_support\n\n                        if (e_support > this.largestEvents_support ) {this.largestEvents_support = e_support;}\n\n                    }\n                }\n\n\n\n            })\n\n\n        }\n\n    }\n\n    get_name(){\n        return this.settings.name\n    }\n\n    set_name(name){\n        this.settings.name = name\n    }\n\n    traverse(o,func_pre, func_post) {\n\n        if (func_pre){\n            func_pre.apply(this,[o,o[\"children\"]])\n        }\n\n        if(o[\"children\"]){\n\n            for (var c in o[\"children\"] ) {\n\n                var child = o[\"children\"][c]\n\n                child = this.traverse(child, func_pre, func_post)\n\n                if (func_post) {\n                    func_post.apply(this,[child,o])\n                }\n\n\n            }\n\n\n        }\n\n        return o\n\n    }\n\n    traverse_hierarchy(o,func_pre, func_post) {\n\n        var children = o[\"children\"] ? o[\"children\"] : o[\"_children\"]\n\n        if (func_pre){\n            func_pre.apply(this,[o,children])\n        }\n\n        if(children ){\n\n            for (var c in children) {\n\n                var child = children[c]\n\n                child = this.traverse_hierarchy(child, func_pre, func_post)\n\n                if (func_post) {\n                    func_post.apply(this,[child,o])\n                }\n\n\n            }\n\n\n        }\n\n        return o\n\n    }\n\n    set_parent(node,parent){\n        node.parent = parent\n    }\n\n    set_cumulated_length(node, children){\n        if (node.parent) {\n            node.depth = node.parent.depth + 1\n            node.distance_to_root = node.parent.distance_to_root + node.branch_length\n        }\n        else{\n            node.distance_to_root = 0\n        node.depth = 0}\n\n    }\n\n    factory(json){ // todo do one traversal with all in one function\n\n        var p;\n\n        //has_branch_lenght\n        this.settings.has_branch_lenght = false;\n\n        json.children.forEach((child) => {\n            if (typeof child.branch_length != 'undefined') { this.settings.has_branch_lenght = true; }\n        })\n\n        if (this.settings.has_branch_lenght) {\n            this.settings.labels['node'].add('Length')\n            this.settings.colorlabels['node'].add('Length')\n            this.settings.extended_data_type['Length'] = 'num'\n\n            this.settings.style.color_extent_max['node']['Length'] = 0;\n            this.settings.style.color_extent_min['node']['Length'] = 1000000000;\n        }\n\n\n        // if branch size is not used put 1\n        if (!this.settings.has_branch_lenght) {\n            p = this.traverse(json, function(n,c){n.branch_length=1})\n            p.branch_length = 0 // root\n        }\n        else{ // sanity check\n            p = this.traverse(json, function(n,c){if (typeof n.branch_length == 'undefined') {n.branch_length=1} })\n            if (typeof p.branch_length == 'undefined') {p.branch_length=1}\n        }\n\n        // set parent attribute\n        p = this.traverse(json, null , this.set_parent)\n\n        // compute cumulated  lenght\n        p = this.traverse(p, this.set_cumulated_length , null)\n\n        this.traverse(p, function(n,c){\n\n            n.extended_informations = {}\n            n.elementS = {}\n            n.elementBCN = {}\n\n            if(n.branch_length){\n                n.extended_informations['Length'] = n.branch_length;\n                if (this.settings.style.color_extent_max['node']['Length'] < n.branch_length){\n                    this.settings.style.color_extent_max['node']['Length'] = n.branch_length\n                }\n\n                if (this.settings.style.color_extent_min['node']['Length'] > n.branch_length){\n                    this.settings.style.color_extent_min['node']['Length'] = n.branch_length\n                }\n            }\n\n            if (typeof c !== 'undefined' && typeof n.name !== 'undefined' && n.name !== \"\" ) {\n                n.extended_informations['Data'] = n.name;\n                this.settings.labels['node'].add('Data')\n                this.settings.extended_data_type['Data'] = 'num'\n\n                if (!isNaN(n.name)){\n\n                    if (!this.settings.colorlabels['node'].has('Data')){\n                        this.settings.colorlabels['node'].add('Data');\n                        this.settings.style.color_extent_max['node']['Data'] = 0;\n                        this.settings.style.color_extent_min['node']['Data'] = 1000000000;\n                    }\n\n                    if (this.settings.style.color_extent_max['node']['Data'] <n.name){\n                        this.settings.style.color_extent_max['node']['Data'] = n.name\n                    }\n\n                    if (this.settings.style.color_extent_min['node']['Data'] > n.name){\n                        this.settings.style.color_extent_min['node']['Data'] = n.name\n                    }\n\n                }\n\n                else {this.settings.extended_data_type['Data'] = 'cat'}\n            }\n\n            if(n.data_nhx && Object.keys(n.data_nhx).length > 0){\n\n                Object.entries(n.data_nhx).forEach(([key, value]) => {\n\n                    this.settings.labels['node'].add(key)\n\n                    switch(key){\n                        case 'Ev':\n                            if (value == 'duplication') {\n                                n.duplication = true\n                                this.settings.has_duplications = true;\n                            }\n                            n.extended_informations.events = value\n                            n.extended_informations[key] = value\n                            this.settings.extended_data_type['Ev'] = 'cat'\n                            break;\n                        case 'DD':\n                        case 'D':\n                            if (value == 'Y') {\n                                n.duplication = true\n                                this.settings.has_duplications = true;\n\n                            }\n                            else if (value == 'N'){\n                                n.duplication = false\n                                this.settings.has_duplications = true;\n                            }\n                            n.extended_informations.events = value\n                            n.extended_informations[key] = value\n                            this.settings.extended_data_type['D'] = 'cat'\n                            break;\n                        default:\n                            n.extended_informations[key] = value\n                            this.settings.extended_data_type[key] = 'cat'\n                            break;\n                    }\n                });\n\n\n            }\n\n\n            if (n.depth > this.settings.tree.max_depth){\n                this.settings.tree.max_depth = n.depth\n            }\n            if (!(n.hasOwnProperty('children'))){\n                this.leaves.push(n)\n                n.correspondingLeaf = {}\n\n            }\n\n            n.leaves = this.get_leaves(n)\n\n        })\n\n        this.settings.suggestions = [] // autocomplete name\n        this.traverse(json, function(n,c){\n            if (n.name !== ''){this.settings.suggestions.push(n.name)}}) //todo add id also and ncBI and more + check empty cfucntion\n        return p\n    }\n\n    build_hierarchy_mockup(){\n        return d3__WEBPACK_IMPORTED_MODULE_0__.hierarchy(this.data, d => d.children );\n    }\n\n    parse(){\n\n        if (this.settings.data_type === \"newick\") {\n            return biojs_io_newick__WEBPACK_IMPORTED_MODULE_2__.parse_newick(this.input_data);\n        }\n\n        else if (this.settings.data_type === \"nhx\") {\n            return parse_nhx(this.input_data);\n        }\n\n        else if (this.settings.data_type === \"json\") {\n            return this.input_data\n        }\n\n\n\n\n    }\n\n    collapse(data, action){\n        if (!data.children){return}\n        if (action) {data.collapse = true}\n        else if (action == false){data.collapse = false}\n        else{data.collapse ? data.collapse = false : data.collapse = true;}\n\n    }\n\n    collapseAll(data, action){\n        if (action) {\n            this.traverse(data, (n,c) => {\n                this.collapse(n, true)\n                //n.collapse = true\n            } , null)}\n        else if (action == false){\n            this.traverse(data, (n,c) => {\n                this.collapse(n, false)\n                //n.collapse = false\n            } , null)}\n\n\n    }\n\n    get_all_collapse(data){\n\n        var collapsed = []\n\n        this.traverse(data, (n,c) => {if (n.collapse){\n                collapsed.push(n)\n            }}, null)\n\n        return collapsed\n\n    }\n\n    apply_collapse_to_list(collapsed){\n        this.traverse(this.data, (n,c) => {\n            if (collapsed.includes(n)){\n                this.collapse(n, true)\n            }\n            else{\n                this.collapse(n, false)\n            }\n        } , null)\n    }\n\n    swap_subtrees(data){\n         var e = data.children.pop()\n        data.children.unshift(e)\n        data.leaves = this.get_leaves(data)\n    }\n\n    unswap_subtrees(data){\n        var e = data.children.shift()\n        data.children.push(e)\n        data.leaves = this.get_leaves(data)\n\n    }\n\n    reroot(data){\n\n        // extract meta data (zoom)\n        var meta = this.zoom;\n\n        // create new root r\n\n        var root = {\"children\": [], \"name\": \"\", \"branch_length\": 0, \"extended_informations\": {}}\n\n        for (var key in data.extended_informations) {\n            if (data.extended_informations.hasOwnProperty(key)) {\n                root.extended_informations[key] =  null\n            }\n        }\n\n        // source and target node of the clicked edges\n        var parent = data.parent\n        var child = data\n\n        // insert new root node between target and source and connect\n        root.children.push(child)\n        parent.children.push(root)\n        this.set_parent(child, root )\n        this.set_parent(root, parent )\n        const index = parent.children.indexOf(child);\n        if (index > -1) {\n            parent.children.splice(index, 1);\n        }\n\n        // ajust distance now that distance target/source is splitted in two\n        var old_distance = child.branch_length\n\n        child.branch_length = old_distance/2\n        child.extended_informations['Length'] = old_distance/2\n\n        parent.branch_length_before_reverse = parent.branch_length\n        parent.branch_length = old_distance /2\n        parent.extended_informations['Length'] = old_distance/2\n\n\n\n        // While we are at the old root reverse child/parent order\n        var parent = parent\n        var child = root\n        var stack = []\n\n        while (parent.root != true) {\n\n            stack.push([parent,child])\n\n            child = parent\n            parent = parent.parent\n\n\n            parent.branch_length_before_reverse = parent.branch_length\n            if (child.branch_length_before_reverse){\n                parent.branch_length = child.branch_length_before_reverse\n                parent.extended_informations['Length'] = child.branch_length_before_reverse\n            }\n            else{\n                parent.branch_length = child.branch_length\n                parent.extended_informations['Length'] = child.branch_length\n            }\n\n\n        }\n        stack.push([parent,child])\n        for (var e in stack){\n            var p = stack[e][0]\n            var c = stack[e][1]\n\n            this.reverse_order(p,c)\n\n        }\n\n        // Remove old root\n\n        var old_root = parent\n        var leading_branch = parent.parent\n\n\n\n        if (old_root.children.length == 1){\n\n            const ce = leading_branch.children.indexOf(old_root);\n            if (ce > -1) {\n                leading_branch.children.splice(ce, 1);\n            }\n\n            var i = 0,len = old_root.children.length;\n            while (i < len) {\n                let c = old_root.children[i]\n                c.parent = leading_branch\n                leading_branch.children.push(c)\n                i++\n            }\n\n            old_root = null\n\n\n\n        }\n\n        // For multifurcation we need to keep the root\n        else {\n            old_root.root = false\n            old_root.branch_length = leading_branch.branch_length\n            parent.extended_informations['Length'] = leading_branch.branch_length\n        }\n\n\n        // configure new root\n        root.zoom = meta\n        this.data = root;\n        this.data.root = true;\n\n        root.leaves = this.get_leaves(root)\n\n\n        this.traverse(root, function(n,c){\n            n.leaves = this.get_leaves(n)\n        })\n\n\n\n    }\n\n    trim(branch){\n\n        // source and target node of the clicked edges\n        var parent = branch.parent\n        var child = branch\n\n        var untrim_data = {\n            \"parent\" : null,\n            \"floating\" : null,\n            \"untrim_data\" : null,\n            \"index\": null,\n            \"root_mode\": false,\n        }\n\n\n        if (parent.children.length > 2) {\n            untrim_data.index = this.detach_child(parent,child)\n\n            untrim_data.parent = parent;\n            untrim_data.floating = false;\n            untrim_data.child =  child;\n\n            return untrim_data;\n        }\n\n        else{\n\n            if (typeof parent.parent == 'undefined'){ // parent is root\n                untrim_data.root_mode = true;\n\n                var sibling = parent.children[0] == child ? parent.children[1] : parent.children[0]\n                untrim_data.index = this.detach_child(parent, sibling)\n\n                this.data = sibling;\n\n                untrim_data.parent = null;\n                untrim_data.floating = parent;\n                untrim_data.child =  sibling;\n\n                return untrim_data;\n\n            }\n\n            else{\n\n                this.detach_child(parent.parent, parent)\n                var sibling = parent.children[0] == child ? parent.children[1] : parent.children[0]\n                untrim_data.index = this.detach_child(parent, sibling)\n                this.attach_child(parent.parent, sibling)\n\n                untrim_data.parent = parent.parent;\n                untrim_data.floating = parent;\n                untrim_data.child =  sibling;\n\n                return untrim_data;\n\n            }\n\n\n        }\n\n    }\n\n    untrim(parent, floating, child, index, root_mode){\n\n        if (floating != false){\n            if (root_mode){\n                this.data = floating\n                this.attach_child(floating,child, index)\n\n            }else{\n                this.detach_child(parent,child)\n                this.attach_child(parent,floating)\n                this.attach_child(floating,child, index)\n            }\n\n        }\n        else {\n            this.attach_child(parent,child, index)\n        }\n    }\n\n    detach_child(parent, child){\n        var index = parent.children.indexOf(child);\n        if (index > -1) {\n            parent.children.splice(index, 1);\n        }\n        return index\n    }\n\n    attach_child(parent,child_to_adopt, index){\n\n        if (typeof index !== 'undefined') {\n            parent.children.splice( index, 0, child_to_adopt );\n            child_to_adopt.parent = parent;\n\n        } else {\n            parent.children.push(child_to_adopt);\n            child_to_adopt.parent = parent;\n        }\n\n\n    }\n\n    interleave_node(parent, to_insert,child){\n        this.detach_child(parent, child)\n        this.attach_child(parent,to_insert)\n        this.attach_child(to_insert, child)\n    }\n\n    store_zoomTransform(zoom){\n\n        this.zoom = {\n            \"k\":zoom.k,\n            \"x\":zoom.x,\n            \"y\":zoom.y,\n        };\n    }\n\n    /**\n     Description:\n     Creates list of leaves of each node in subtree rooted at v\n\n     Note:\n     Difference between deep leaf list and leaves in:\n     (A:0.1,B:0.2,(C:0.3,D:0.4):0.5);\n     - Root has leaves: A, B, C and D (terminal leaves)\n     - Root has deep leaves: A, B, C, D and CD (terminal leaves + intermediate leaves)\n     */\n    createDeepLeafList(filter) {\n\n        function is_leaf(str) {\n            return !str.includes(\"||\");\n        }\n\n         var build_deepLeafList = function(child, node){\n\n             if ( child.hasOwnProperty('children') ){\n                 var dp = child.deepLeafList.filter(is_leaf).sort()\n                 if (!dp.every((e) => e === '')){\n                     child.deepLeafList.push(dp.join('||'));\n                 }\n\n             }\n\n             node.deepLeafList = node.deepLeafList.concat(child.deepLeafList)\n        }\n\n        var build_deepLeafLeaves = function(node,children){\n\n             if (!(node.hasOwnProperty('children') )){\n\n                 if (typeof filter != 'undefined') {\n                     if (filter.includes(node.name)){\n                         node.deepLeafList = [node.name]\n                     }\n                     else{\n                         node.deepLeafList = []\n                     }\n\n                 }\n                 else{\n                     node.deepLeafList = [node.name]\n                 }\n\n             }\n             else {\n                 node.deepLeafList = []\n             }\n\n\n        }\n\n        this.traverse(this.data, build_deepLeafLeaves, build_deepLeafList)\n\n    }\n\n    createMinHash(){\n\n        var assign_hash = function(node,children){\n\n            node.min_hash = new minhashjs__WEBPACK_IMPORTED_MODULE_1__.MinHash.MinHash()\n            node.deepLeafList.map(function(w) { node.min_hash.update(w) });\n        }\n\n        this.traverse(this.data, assign_hash, null)\n    }\n\n    removeMinHash(){\n\n        var remove_hash = function(node,children){\n\n            node.min_hash = null\n        }\n\n        this.traverse(this.data, remove_hash, null)\n\n    }\n\n    reverse_order(parent,child) {\n\n        child.children.push(parent)\n        parent.parent =child\n\n        const b = parent.children.indexOf(child);\n        if (b > -1) {\n            parent.children.splice(b, 1);\n        }\n\n    }\n\n    get_leaves(node){\n\n\n        var l = []\n\n        this.traverse(node, function(n,c){\n            if (!(n.hasOwnProperty('children'))){\n                l.push(n)\n            }\n\n\n    })\n        return l\n    }\n\n    remove_circularity(){ // safe my model\n        var data = Object.assign({}, this.data);\n\n        this.traverse(data, function(n,c){\n            n.parent=null;\n            n.leaves=null;\n            n.correspondingLeaf = {}\n            n.elementBCN = null})\n\n        return data\n    }\n\n    remove_circularity_only_parent_and_leaves(){ // safe my model\n        var data = Object.assign({}, this.data);\n\n        this.traverse(data, function(n,c){\n            n.parent=null;\n            n.leaves=null;\n        })\n\n        return data\n    }\n\n    add_circularity_back(){\n\n        this.data.leaves = this.get_leaves(this.data)\n\n        this.traverse(this.data, function(n,c){n.leaves = this.get_leaves(n)}, this.set_parent)\n\n\n    }\n\n    add_meta_leaves(meta, headers, api){\n\n        // headers: column_name -> type\n\n        Object.keys(headers).forEach(item => {\n            if (item != 'id' || item != 'Length' ) {\n\n                this.settings.extended_data_type[item] = headers[item]\n                this.settings.domain_extended_data[item] = []\n                this.settings.labels['leaf'].add(item)\n                this.settings.colorlabels['leaf'].add(item)\n\n                if (headers[item] == 'num'){\n                    this.settings.style.color_extent_max['leaf'][item] = 0\n                    this.settings.style.color_extent_min['leaf'][item] = 100000\n                }\n            }\n\n        })\n\n            this.get_leaves(this.data).forEach(d => {\n            if (d.name in meta){\n\n                Object.entries(meta[d.name]).forEach(item => {\n                    if (item[0] != 'id'){\n\n                        d.extended_informations[item[0]]= item[1]\n\n                        if (this.settings.extended_data_type[item[0]] == 'num' && !isNaN(item[1])) {\n\n                            item[1] = item[1].toString().indexOf('.') != -1 ? parseFloat(item[1]) : parseInt(item[1])\n\n                            if (this.settings.style.color_extent_max['leaf'][item[0]] < item[1]) {\n                                this.settings.style.color_extent_max['leaf'][item[0]] = item[1]\n                            }\n\n                            if (this.settings.style.color_extent_min['leaf'][item[0]] > item[1]) {\n                                this.settings.style.color_extent_min['leaf'][item[0]] = item[1]\n                            }\n\n                        }\n\n                        if (this.settings.extended_data_type[item[0]] == 'cat'){\n\n                            var cs = api.get_color_scale(item[0])\n                            cs.add_value_to_map(item[1])\n\n                            this.settings.domain_extended_data[item[0]].push(item[1])\n                        }\n\n                    }\n\n                })\n\n            }\n\n        })\n\n\n        Object.keys(headers).forEach(item => {\n            if (item != 'id' || item != 'Length' ) {\n\n\n                if (headers[item] == 'cat'){\n                    api.get_color_scale(item).update()\n\n                }\n            }\n\n        })\n\n\n\n    }\n\n    add_meta_nodes(meta, headers, api){\n\n\n        Object.keys(headers).forEach(item => {\n            if (item != 'id' || item != 'Length' ) {\n                this.settings.extended_data_type[item] = headers[item]\n                this.settings.domain_extended_data[item] = []\n                this.settings.labels['node'].add(item)\n                this.settings.colorlabels['node'].add(item)\n\n                if (headers[item] == 'num'){\n                    this.settings.style.color_extent_max['node'][item] = 0\n                    this.settings.style.color_extent_min['node'][item] = 100000\n                }\n            }\n\n        })\n\n\n        this.traverse(this.data, function(n,c){\n\n            if (n.extended_informations['Data'] in meta){\n\n                Object.entries(meta[n.extended_informations['Data']]).forEach(item => {\n                    if (item[0] != 'id'){\n                        n.extended_informations[item[0]]= item[1]\n\n                        if (this.settings.extended_data_type[item[0]] == 'cat'){\n\n                            var cs = api.get_color_scale(item[0])\n                            cs.add_value_to_map(item[1])\n\n                        }\n\n                        else if (this.settings.extended_data_type[item[0]] == 'num' && !isNaN(item[1]) ) {\n\n                            item[1] = item[1].toString().indexOf('.') != -1 ? parseFloat(item[1]) : parseInt(item[1])\n\n                            if (this.settings.style.color_extent_max['node'][item[0]] < item[1]) {\n                                this.settings.style.color_extent_max['node'][item[0]] = item[1]\n                            }\n\n                            if (this.settings.style.color_extent_min['node'][item[0]] > item[1]) {\n                                this.settings.style.color_extent_min['node'][item[0]] = item[1]\n                            }\n\n                        }\n                    }\n                })\n\n            }\n        })\n\n        Object.keys(headers).forEach(item => {\n\n            if (item != 'id' || item != 'Length' ) {\n\n\n                if (headers[item] == 'cat'){\n                    api.get_color_scale(item).update()\n\n                }\n            }\n\n        })\n    }\n\n    get_node_by_leafset(lset){\n\n        function setsAreEqual(a, b) {\n            if (a.size !== b.size) {\n                return false;\n            }\n\n            return Array.from(a).every(element => {\n                return b.has(element);\n            });\n        }\n\n        lset = new Set(lset.map(leaf => leaf.toString()))\n\n        var target = false\n\n\n        var check = function(node,children){\n\n            var nl = new Set(node.leaves.map(leaf => leaf.name.replaceAll(\"'\", '').toString()))\n\n            if ( setsAreEqual(nl,lset)){\n                target = node\n            }\n\n        }\n\n        this.traverse(this.data, check, null)\n\n        return target\n    }\n\n};\n\n//# sourceURL=webpack://PhyloIO/./src/model.js?");

/***/ }),

/***/ "./src/worker_bcn.js":
/*!***************************!*\
  !*** ./src/worker_bcn.js ***!
  \***************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
eval("__webpack_require__.r(__webpack_exports__);\n/* harmony import */ var minhashjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! minhashjs */ \"./node_modules/minhashjs/index.js\");\n/* harmony import */ var _model__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./model */ \"./src/model.js\");\n\n\n\n\nself.onmessage = (event) => {\n\n    var containers = event.data;\n\n    var processed_tree = compute_similarity_container_pair(containers.tree1,containers.tree2)\n\n    postMessage(processed_tree);\n\n    self.close();\n};\n\n\n\nfunction compute_similarity_container_pair(t1,t2){\n\n    var t1 = new _model__WEBPACK_IMPORTED_MODULE_1__[\"default\"](t1.data, t1.settings, false)\n    var t2 = new _model__WEBPACK_IMPORTED_MODULE_1__[\"default\"](t2.data, t2.settings, false)\n\n\n    console.time(\"similarity\");\n\n    // X = Intersection of T1 & T2 leaves.\n    if (t1.leaves.length <= 0){\n        t1.leaves = t1.get_leaves(t1.data)\n    }\n    if (t2.leaves.length <= 0){\n        t2.leaves = t2.get_leaves(t2.data)\n    }\n    var common_leaves = t1.leaves.map(leaf => leaf.name).filter(value => t2.leaves.map(leaf => leaf.name).includes(value));\n    console.log(\"Intersection\");\n    console.timeLog(\"similarity\");\n\n    // DeepLeaf with filter(X = True)\n    t1.createDeepLeafList(common_leaves)\n    t2.createDeepLeafList(common_leaves)\n    console.log(\"DeepLeaf\");\n    console.timeLog(\"similarity\");\n\n    // MinHash\n    t1.createMinHash()\n    t2.createMinHash()\n    console.log(\"MinHash\");\n    console.timeLog(\"similarity\");\n\n    // For all T1 & T2 Nodes compute the BCN with MinHash\n    var nodes_t1 = []\n    var nodes_t2 = []\n\n    var forest1 = new minhashjs__WEBPACK_IMPORTED_MODULE_0__.MinHashLSHForest.MinHashLSHForest()\n    var forest2 = new minhashjs__WEBPACK_IMPORTED_MODULE_0__.MinHashLSHForest.MinHashLSHForest()\n\n    var cpt =0\n    t1.traverse(t1.data, function(h,children){nodes_t1.push(h);forest1.add(h, h.min_hash);cpt++}, null)\n    t2.traverse(t2.data, function(z,children){nodes_t2.push(z);forest2.add(z, z.min_hash);cpt++}, null)\n\n    forest1.index()\n    forest2.index()\n\n    find_BCN(nodes_t1, forest2, t2.uid)\n    find_BCN(nodes_t2, forest1, t1.uid)\n\n    console.log(\"BCN done\");\n\n    console.timeLog(\"similarity\");\n\n    // Clean non essential datume\n    t1.removeMinHash()\n    t2.removeMinHash()\n\n    console.timeEnd(\"similarity\");\n\n    t1.settings.similarity.push(t2.uid)\n    t2.settings.similarity.push(t1.uid)\n\n    t1.data = t1.remove_circularity_only_parent_and_leaves()\n    t2.data = t2.remove_circularity_only_parent_and_leaves()\n\n    return [t1,t2]\n\n}\n\nfunction find_BCN(nodes_list, target_forest, target_uid){\n    nodes_list.forEach((node) => {\n\n        function is_leaf(str) {\n            return !str.includes(\"||\");\n        }\n\n        var matches = target_forest.query(node.min_hash,10)\n\n        var l = new Set(node.deepLeafList.filter(is_leaf))\n        if (l.size > 0){\n\n            var max_jacc = 0\n            var BCN = null\n\n            matches.forEach(e => {\n                var r =   new Set(e.deepLeafList.filter(is_leaf))\n\n                if (r.size > 0){\n\n\n                    var inter = Array.from(r).filter(x => l.has(x)).length\n                    var union = [...new Set([...l, ...r])].length;\n\n                    var jj = inter/union\n\n\n\n                    if (jj > max_jacc){\n                        max_jacc = jj\n                        BCN = e\n                    }\n\n                }\n\n\n\n\n\n            })\n\n            if (max_jacc > 0) {\n                node.elementS[target_uid] = max_jacc\n                if (!node.elementBCN){\n                    node.elementBCN = {}\n                }\n                node.elementBCN[target_uid] = BCN\n            }\n\n\n        }\n\n\n\n\n    })\n}\n\n//# sourceURL=webpack://PhyloIO/./src/worker_bcn.js?");

/***/ }),

/***/ "?d546":
/*!************************!*\
  !*** buffer (ignored) ***!
  \************************/
/***/ (() => {

eval("/* (ignored) */\n\n//# sourceURL=webpack://PhyloIO/buffer_(ignored)?");

/***/ }),

/***/ "?8131":
/*!************************!*\
  !*** buffer (ignored) ***!
  \************************/
/***/ (() => {

eval("/* (ignored) */\n\n//# sourceURL=webpack://PhyloIO/buffer_(ignored)?");

/***/ }),

/***/ "?3fc0":
/*!************************!*\
  !*** crypto (ignored) ***!
  \************************/
/***/ (() => {

eval("/* (ignored) */\n\n//# sourceURL=webpack://PhyloIO/crypto_(ignored)?");

/***/ }),

/***/ "?4068":
/*!************************!*\
  !*** buffer (ignored) ***!
  \************************/
/***/ (() => {

eval("/* (ignored) */\n\n//# sourceURL=webpack://PhyloIO/buffer_(ignored)?");

/***/ }),

/***/ "?e7e4":
/*!************************!*\
  !*** buffer (ignored) ***!
  \************************/
/***/ (() => {

eval("/* (ignored) */\n\n//# sourceURL=webpack://PhyloIO/buffer_(ignored)?");

/***/ }),

/***/ "?7bec":
/*!************************!*\
  !*** buffer (ignored) ***!
  \************************/
/***/ (() => {

eval("/* (ignored) */\n\n//# sourceURL=webpack://PhyloIO/buffer_(ignored)?");

/***/ }),

/***/ "?0aec":
/*!************************!*\
  !*** buffer (ignored) ***!
  \************************/
/***/ (() => {

eval("/* (ignored) */\n\n//# sourceURL=webpack://PhyloIO/buffer_(ignored)?");

/***/ }),

/***/ "?fbf1":
/*!************************!*\
  !*** buffer (ignored) ***!
  \************************/
/***/ (() => {

eval("/* (ignored) */\n\n//# sourceURL=webpack://PhyloIO/buffer_(ignored)?");

/***/ }),

/***/ "?ed1b":
/*!**********************!*\
  !*** util (ignored) ***!
  \**********************/
/***/ (() => {

eval("/* (ignored) */\n\n//# sourceURL=webpack://PhyloIO/util_(ignored)?");

/***/ }),

/***/ "?d17e":
/*!**********************!*\
  !*** util (ignored) ***!
  \**********************/
/***/ (() => {

eval("/* (ignored) */\n\n//# sourceURL=webpack://PhyloIO/util_(ignored)?");

/***/ })

/******/ 	});
/************************************************************************/
/******/ 	// The module cache
/******/ 	var __webpack_module_cache__ = {};
/******/ 	
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/ 		// Check if module is in cache
/******/ 		var cachedModule = __webpack_module_cache__[moduleId];
/******/ 		if (cachedModule !== undefined) {
/******/ 			return cachedModule.exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = __webpack_module_cache__[moduleId] = {
/******/ 			id: moduleId,
/******/ 			loaded: false,
/******/ 			exports: {}
/******/ 		};
/******/ 	
/******/ 		// Execute the module function
/******/ 		__webpack_modules__[moduleId].call(module.exports, module, module.exports, __webpack_require__);
/******/ 	
/******/ 		// Flag the module as loaded
/******/ 		module.loaded = true;
/******/ 	
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/ 	
/******/ 	// expose the modules object (__webpack_modules__)
/******/ 	__webpack_require__.m = __webpack_modules__;
/******/ 	
/******/ 	// the startup function
/******/ 	__webpack_require__.x = () => {
/******/ 		// Load entry module and return exports
/******/ 		// This entry module depends on other loaded chunks and execution need to be delayed
/******/ 		var __webpack_exports__ = __webpack_require__.O(undefined, ["vendors-node_modules_d3_index_js-node_modules_file-saver_dist_FileSaver_min_js","vendors-node_modules_biojs-io-newick_src_index_js-node_modules_minhashjs_index_js","src_utils_js"], () => (__webpack_require__("./src/worker_bcn.js")))
/******/ 		__webpack_exports__ = __webpack_require__.O(__webpack_exports__);
/******/ 		return __webpack_exports__;
/******/ 	};
/******/ 	
/************************************************************************/
/******/ 	/* webpack/runtime/chunk loaded */
/******/ 	(() => {
/******/ 		var deferred = [];
/******/ 		__webpack_require__.O = (result, chunkIds, fn, priority) => {
/******/ 			if(chunkIds) {
/******/ 				priority = priority || 0;
/******/ 				for(var i = deferred.length; i > 0 && deferred[i - 1][2] > priority; i--) deferred[i] = deferred[i - 1];
/******/ 				deferred[i] = [chunkIds, fn, priority];
/******/ 				return;
/******/ 			}
/******/ 			var notFulfilled = Infinity;
/******/ 			for (var i = 0; i < deferred.length; i++) {
/******/ 				var [chunkIds, fn, priority] = deferred[i];
/******/ 				var fulfilled = true;
/******/ 				for (var j = 0; j < chunkIds.length; j++) {
/******/ 					if ((priority & 1 === 0 || notFulfilled >= priority) && Object.keys(__webpack_require__.O).every((key) => (__webpack_require__.O[key](chunkIds[j])))) {
/******/ 						chunkIds.splice(j--, 1);
/******/ 					} else {
/******/ 						fulfilled = false;
/******/ 						if(priority < notFulfilled) notFulfilled = priority;
/******/ 					}
/******/ 				}
/******/ 				if(fulfilled) {
/******/ 					deferred.splice(i--, 1)
/******/ 					var r = fn();
/******/ 					if (r !== undefined) result = r;
/******/ 				}
/******/ 			}
/******/ 			return result;
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/define property getters */
/******/ 	(() => {
/******/ 		// define getter functions for harmony exports
/******/ 		__webpack_require__.d = (exports, definition) => {
/******/ 			for(var key in definition) {
/******/ 				if(__webpack_require__.o(definition, key) && !__webpack_require__.o(exports, key)) {
/******/ 					Object.defineProperty(exports, key, { enumerable: true, get: definition[key] });
/******/ 				}
/******/ 			}
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/ensure chunk */
/******/ 	(() => {
/******/ 		__webpack_require__.f = {};
/******/ 		// This file contains only the entry chunk.
/******/ 		// The chunk loading function for additional chunks
/******/ 		__webpack_require__.e = (chunkId) => {
/******/ 			return Promise.all(Object.keys(__webpack_require__.f).reduce((promises, key) => {
/******/ 				__webpack_require__.f[key](chunkId, promises);
/******/ 				return promises;
/******/ 			}, []));
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/get javascript chunk filename */
/******/ 	(() => {
/******/ 		// This function allow to reference async chunks and sibling chunks for the entrypoint
/******/ 		__webpack_require__.u = (chunkId) => {
/******/ 			// return url for filenames based on template
/******/ 			return "" + chunkId + ".phylo.js";
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/global */
/******/ 	(() => {
/******/ 		__webpack_require__.g = (function() {
/******/ 			if (typeof globalThis === 'object') return globalThis;
/******/ 			try {
/******/ 				return this || new Function('return this')();
/******/ 			} catch (e) {
/******/ 				if (typeof window === 'object') return window;
/******/ 			}
/******/ 		})();
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/hasOwnProperty shorthand */
/******/ 	(() => {
/******/ 		__webpack_require__.o = (obj, prop) => (Object.prototype.hasOwnProperty.call(obj, prop))
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/make namespace object */
/******/ 	(() => {
/******/ 		// define __esModule on exports
/******/ 		__webpack_require__.r = (exports) => {
/******/ 			if(typeof Symbol !== 'undefined' && Symbol.toStringTag) {
/******/ 				Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
/******/ 			}
/******/ 			Object.defineProperty(exports, '__esModule', { value: true });
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/node module decorator */
/******/ 	(() => {
/******/ 		__webpack_require__.nmd = (module) => {
/******/ 			module.paths = [];
/******/ 			if (!module.children) module.children = [];
/******/ 			return module;
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/publicPath */
/******/ 	(() => {
/******/ 		__webpack_require__.p = "";
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/importScripts chunk loading */
/******/ 	(() => {
/******/ 		// no baseURI
/******/ 		
/******/ 		// object to store loaded chunks
/******/ 		// "1" means "already loaded"
/******/ 		var installedChunks = {
/******/ 			"src_worker_bcn_js": 1
/******/ 		};
/******/ 		
/******/ 		// importScripts chunk loading
/******/ 		var installChunk = (data) => {
/******/ 			var [chunkIds, moreModules, runtime] = data;
/******/ 			for(var moduleId in moreModules) {
/******/ 				if(__webpack_require__.o(moreModules, moduleId)) {
/******/ 					__webpack_require__.m[moduleId] = moreModules[moduleId];
/******/ 				}
/******/ 			}
/******/ 			if(runtime) runtime(__webpack_require__);
/******/ 			while(chunkIds.length)
/******/ 				installedChunks[chunkIds.pop()] = 1;
/******/ 			parentChunkLoadingFunction(data);
/******/ 		};
/******/ 		__webpack_require__.f.i = (chunkId, promises) => {
/******/ 			// "1" is the signal for "already loaded"
/******/ 			if(!installedChunks[chunkId]) {
/******/ 				if(true) { // all chunks have JS
/******/ 					importScripts(__webpack_require__.p + __webpack_require__.u(chunkId));
/******/ 				}
/******/ 			}
/******/ 		};
/******/ 		
/******/ 		var chunkLoadingGlobal = self["webpackChunkPhyloIO"] = self["webpackChunkPhyloIO"] || [];
/******/ 		var parentChunkLoadingFunction = chunkLoadingGlobal.push.bind(chunkLoadingGlobal);
/******/ 		chunkLoadingGlobal.push = installChunk;
/******/ 		
/******/ 		// no HMR
/******/ 		
/******/ 		// no HMR manifest
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/startup chunk dependencies */
/******/ 	(() => {
/******/ 		var next = __webpack_require__.x;
/******/ 		__webpack_require__.x = () => {
/******/ 			return Promise.all(["vendors-node_modules_d3_index_js-node_modules_file-saver_dist_FileSaver_min_js","vendors-node_modules_biojs-io-newick_src_index_js-node_modules_minhashjs_index_js","src_utils_js"].map(__webpack_require__.e, __webpack_require__)).then(next);
/******/ 		};
/******/ 	})();
/******/ 	
/************************************************************************/
/******/ 	
/******/ 	// run startup
/******/ 	var __webpack_exports__ = __webpack_require__.x();
/******/ 	PhyloIO = __webpack_exports__;
/******/ 	
/******/ })()
;