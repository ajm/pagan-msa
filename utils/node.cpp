#include "node.h"
#include <iostream>

using namespace std;
using namespace ppa;

//Node::Node()
//{
//    isleaf = true;
//}

Node::~Node()
{
    if(sequence!=0)
        delete sequence;
}

void Node::add_sequence( string seq_string, string full_dna_alphabet)
{
    if(Settings::noise>4)
        cout<<"Node::add_sequence "<<name<<"\n";
    sequence = new Sequence(seq_string, full_dna_alphabet);
}

void Node::get_alignment(vector<Fasta_entry> *aligned_sequences,bool include_internal_nodes)
{
    vector<Node*> nodes;
    if(include_internal_nodes)
        this->get_all_nodes(&nodes);
    else
        this->get_leaf_nodes(&nodes);

    for(unsigned int i=0;i<nodes.size();i++)
    {
        Fasta_entry entry;
        entry.name = nodes.at(i)->get_name();
        entry.comment = nodes.at(i)->get_name_comment();

        aligned_sequences->push_back(entry);
    }

    Sequence *root = this->get_sequence();
    int root_length = root->sites_length();

    for(int j=1;j<root_length-1;j++)
    {
        vector<char> column;
        this->get_alignment_column_at(j,&column,include_internal_nodes);

        for(unsigned int i=0;i<aligned_sequences->size();i++)
        {
            aligned_sequences->at(i).sequence.push_back(column.at(i));
        }
    }
}

void Node::get_alignment_column_at(int j,vector<char> *column,bool include_internal_nodes)
{
    if(leaf)
    {
        int state = sequence->get_site_at(j)->get_state();
        column->push_back(sequence->get_full_alphabet().at(state));
    }
    else
    {
        Site_children *offspring = sequence->get_site_at(j)->get_children();
        int lj = offspring->left_index;
        if(lj>=0)
        {
            left_child->get_alignment_column_at(lj,column,include_internal_nodes);
        }
        else
        {
            int nl = left_child->get_number_of_leaves();
            if(include_internal_nodes)
                nl = left_child->get_number_of_nodes();

            for(int i=0;i<nl;i++)
                column->push_back('-');
        }

        if(include_internal_nodes)
        {
            int cstate = sequence->get_site_at(j)->get_state();
            char c = sequence->get_full_alphabet().at(cstate);

            int pstate = sequence->get_site_at(j)->get_path_state();
            if( pstate == Site::xskipped || pstate == Site::yskipped )
                c = '-';
            column->push_back(c);
        }

        int rj = offspring->right_index;
        if(rj>=0)
        {
            right_child->get_alignment_column_at(rj,column,include_internal_nodes);
        }
        else
        {
            int nl = right_child->get_number_of_leaves();
            if(include_internal_nodes)
                nl = right_child->get_number_of_nodes();

            for(int i=0;i<nl;i++)
                column->push_back('-');
        }
    }
}


void Node::write_metapost_sequence_graph(ostream *output, ostream *output2, int *count, int root_length) const throw (Exception)
{
    *output<<"beginfig("<<*count<<");\npickup pencircle scaled 1pt;\npath c[];\ndefaultscale := 0.5;\n";
    vector<Site> *sites = this->sequence->get_sites();
    string full_alphabet = this->sequence->get_full_alphabet();

    stringstream all_chars;
    for(unsigned int i=0;i<sites->size();i++)
    {
        Site *tsite =  &sites->at(i);

        char c = 's';
        if(tsite->get_site_type()==Site::real_site)
            c = full_alphabet.at(tsite->get_state());
        else if(tsite->get_site_type()==Site::stop_site)
            c = 'e';

        string color = this->get_node_fill_color(c);
        if(tsite->get_branch_count_since_last_used()>0)
            color = "0.5white";

        *output<<"c"<<i<<" = circle"<<"(("<<0.5*i<<"cm,0cm),\""<<c<<"\","<<color<<");\n";
    }

    if(leaf)
        *output<<"label.top(btex $"<<this->get_name()<<"$ etex,(0.125cm,0.25cm));\n";
    else
    {
        string n = this->get_name();
        n = n.substr(1,n.length()-2);
        *output<<"label.top(btex \\#"<<n<<"\\# etex,(0.125cm,0.25cm));\n";
    }

    *output<<"defaultscale := 0.25;\n";

    for(unsigned int i=1;i<sites->size();i++)
    {
        Site *tsite =  &sites->at(i);

        if(tsite->has_bwd_edge())
        {
            Edge *tedge = tsite->get_first_bwd_edge();
            int start = tedge->get_start_site_index();
            int stop  = tedge->get_end_site_index();

            int angle = 0;
            string place = "edgetop";
            if(start+1==stop)
                place = "edgebot";
            else if(start+2==stop)
                angle = 40;
            else if(start+3==stop)
                angle = 30;
            else if(start+4<=stop)
                angle = 20;

            stringstream label;
            if(tedge->get_branch_count_since_last_used()>0)
                label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

            *output<<place<<"(c"<<start<<",c"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            while(tsite->has_next_bwd_edge())
            {
                tedge = tsite->get_next_bwd_edge();
                start = tedge->get_start_site_index();
                stop  = tedge->get_end_site_index();

                angle = 0;
                place = "edgetop";
                if(start+1==stop)
                    place = "edgebot";
                else if(start+2==stop)
                    angle = 40;
                else if(start+3==stop)
                    angle = 30;
                else if(start+4<=stop)
                    angle = 20;

                label.str("");

                if(tedge->get_branch_count_since_last_used()>0)
                    label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

                *output<<place<<"(c"<<start<<",c"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            }
        }
    }

    *output<<"endfig;\n";


    string file = Settings_handle::st.get("mpost-graph-file").as<string>();
    float width = (float)this->sequence->get_sites()->size()/(float)root_length;

    *output2<<"\\includegraphics[width="<<width<<"\\columnwidth]{"<<file<<"."<<*count<<"}\n\n\\bigskip\n";

    ++*count;

}

void Node::write_metapost_alignment_graph(ostream *output, ostream *output2, int *count, int root_length) const throw (Exception)
{
    vector<Site> *sites = this->sequence->get_sites();

    vector<int> left_child_index;
    vector<int> right_child_index;

    for(unsigned int i=0;i<sites->size();i++)
    {
        Site_children *offspring = sites->at(i).get_children();

        if(offspring->left_index>=0)
        {
            left_child_index.push_back(i);
        }
        if(offspring->right_index>=0)
        {
            right_child_index.push_back(i);
        }
    }


    *output<<"beginfig("<<*count<<");\npickup pencircle scaled 1pt;\npath l[]; path r[];\ndefaultscale := 0.5;\n";
    string full_alphabet = this->sequence->get_full_alphabet();

    *output<<"l0 = circle((0cm,1.5cm),\"s\",white);\n";
    *output<<"r0 = circle((0cm,0cm),\"s\",white);\n";


    stringstream all_chars;
    for(unsigned int i=1;i<sites->size();i++)
    {
        Site_children *offspring = sites->at(i).get_children();

        if(offspring->left_index>=0)
        {
            Site *lsite = this->left_child->get_sequence()->get_site_at(offspring->left_index);

            char lc = 's';
            if(lsite->get_site_type()==Site::real_site)
                lc = full_alphabet.at(lsite->get_state());
            else if(lsite->get_site_type()==Site::stop_site)
                lc = 'e';

            string color = this->get_node_fill_color(lc);
            if(lsite->get_branch_count_since_last_used()>0)
                color = "0.5white";

            *output<<"l"<<i<<" = circle"<<"(("<<0.5*i<<"cm,1.5cm),\""<<lc<<"\","<<color<<");\n";
        }
        if(offspring->right_index>=0)
        {
            Site *rsite = this->right_child->get_sequence()->get_site_at(offspring->right_index);

            char rc = 's';
            if(rsite->get_site_type()==Site::real_site)
                rc = full_alphabet.at(rsite->get_state());
            else if(rsite->get_site_type()==Site::stop_site)
                rc = 'e';

            string color = this->get_node_fill_color(rc);
            if(rsite->get_branch_count_since_last_used()>0)
                color = "0.5white";

            *output<<"r"<<i<<" = circle"<<"(("<<0.5*i<<"cm,0cm),\""<<rc<<"\","<<color<<");\n";
        }
    }


    if(left_child->is_leaf())
        *output<<"label.top(btex $"<<left_child->get_name()<<"$ etex,(0.125cm,1.75cm));\n";
    else
    {
        string n = left_child->get_name();
        n = n.substr(1,n.length()-2);
        *output<<"label.top(btex \\#"<<n<<"\\# etex,(0.125cm,1.75cm));\n";
    }

    if(right_child->is_leaf())
        *output<<"label.top(btex $"<<right_child->get_name()<<"$ etex,(0.125cm,0.25cm));\n";
    else
    {
        string n = right_child->get_name();
        n = n.substr(1,n.length()-2);
        *output<<"label.top(btex \\#"<<n<<"\\# etex,(0.125cm,0.25cm));\n";
    }

    *output<<"defaultscale := 0.25;\n";

    for(unsigned int i=1;i<sites->size();i++)
    {

        Site_children *offspring = sites->at(i).get_children();

        if(offspring->left_index>=0)
        {
            Site *tsite = left_child->get_sequence()->get_site_at(offspring->left_index);

            if(tsite->has_bwd_edge())
            {
                Edge *tedge = tsite->get_first_bwd_edge();
                int start = tedge->get_start_site_index();
                int stop  = tedge->get_end_site_index();

                int angle = 0;
                string place = "edgetop";
                if(start+1==stop)
                    place = "edgebot";
                else if(start+2==stop)
                    angle = 40;
                else if(start+3==stop)
                    angle = 30;
                else if(start+4<=stop)
                    angle = 20;

                start = left_child_index.at( start );
                stop  = left_child_index.at( stop );

                stringstream label;
                if(tedge->get_branch_count_since_last_used()>0)
                    label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

                *output<<place<<"(l"<<start<<",l"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

                while(tsite->has_next_bwd_edge())
                {
                    tedge = tsite->get_next_bwd_edge();
                    start = tedge->get_start_site_index();
                    stop  = tedge->get_end_site_index();

                    angle = 0;
                    place = "edgetop";
                    if(start+1==stop)
                        place = "edgebot";
                    else if(start+2==stop)
                        angle = 40;
                    else if(start+3==stop)
                        angle = 30;
                    else if(start+4<=stop)
                        angle = 20;

                    start = left_child_index.at( start );
                    stop  = left_child_index.at( stop );

                    label.str("");

                    if(tedge->get_branch_count_since_last_used()>0)
                        label<<"["<<tedge->get_branch_count_since_last_used()<<"~"<<tedge->get_branch_count_as_skipped_edge()<<"~"<<tedge->get_branch_distance_since_last_used()<<"]";

                    *output<<place<<"(l"<<start<<",l"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

                }
            }
        }

        if(offspring->right_index>=0)
        {
            Site *tsite = right_child->get_sequence()->get_site_at(offspring->right_index);

            if(tsite->has_bwd_edge())
            {
                Edge *tedge = tsite->get_first_bwd_edge();
                int start = tedge->get_start_site_index();
                int stop  = tedge->get_end_site_index();

//                int angle = 0;
//                string place = "edgebot";
//                if(start+1==stop)
//                    place = "edgetop";
//                else if(start+2==stop)
//                    angle = 320;
//                else if(start+3==stop)
//                    angle = 330;
//                else if(start+4<=stop)
//                    angle = 340;
                int angle = 0;
                string place = "edgetop";
                if(start+1==stop)
                    place = "edgebot";
                else if(start+2==stop)
                    angle = 40;
                else if(start+3==stop)
                    angle = 30;
                else if(start+4<=stop)
                    angle = 20;

                start = right_child_index.at( start );
                stop  = right_child_index.at( stop );

                stringstream label;
                if(tedge->get_branch_count_since_last_used()>0)
                    label<<"["<<tedge->get_branch_count_since_last_used()<<"~"<<tedge->get_branch_count_as_skipped_edge()<<"~"<<tedge->get_branch_distance_since_last_used()<<"]";

                *output<<place<<"(r"<<start<<",r"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

                while(tsite->has_next_bwd_edge())
                {
                    tedge = tsite->get_next_bwd_edge();
                    start = tedge->get_start_site_index();
                    stop  = tedge->get_end_site_index();

//                    angle = 0;
//                    place = "edgebot";
//                    if(start+1==stop)
//                        place = "edgetop";
//                    else if(start+2==stop)
//                        angle = 320;
//                    else if(start+3==stop)
//                        angle = 330;
//                    else if(start+4<=stop)
//                        angle = 340;
                    angle = 0;
                    place = "edgetop";
                    if(start+1==stop)
                        place = "edgebot";
                    else if(start+2==stop)
                        angle = 40;
                    else if(start+3==stop)
                        angle = 30;
                    else if(start+4<=stop)
                        angle = 20;

                    start = right_child_index.at( start );
                    stop  = right_child_index.at( stop );

                    label.str("");

                    if(tedge->get_branch_count_since_last_used()>0)
                        label<<"["<<tedge->get_branch_count_since_last_used()<<"~"<<tedge->get_branch_count_as_skipped_edge()<<"~"<<tedge->get_branch_distance_since_last_used()<<"]";

                    *output<<place<<"(r"<<start<<",r"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

                }
            }
        }
    }

    *output<<"endfig;\n";

    string file = Settings_handle::st.get("mpost-graph-file").as<string>();
    float width = (float)this->sequence->get_sites()->size()/(float)root_length;

    *output2<<"\\includegraphics[width="<<width<<"\\columnwidth]{"<<file<<"."<<*count<<"}\n\n\\bigskip\n";
    *output2<<"~\n\n\\bigskip\n";

    ++*count;
}

void Node::check_valid_graph() const
{
    vector<Site> *sites = sequence->get_sites();

    for(unsigned int i=0;i<sites->size();i++)
    {
        Site *ssite = &sites->at(i);
        if( ssite->has_fwd_edge() )
        {
            Edge *edge = ssite->get_first_fwd_edge();
            Site *esite = &sites->at(edge->get_end_site_index());

            if(!esite->contains_bwd_edge(edge,true))
            {
                cout<<"site "<<i<<" has fwd edge from "<<edge->get_start_site_index()<<" to "
                        <<edge->get_end_site_index()<<" but no return\n";
            }

            while( ssite->has_next_fwd_edge() )
            {
                edge = ssite->get_next_fwd_edge();
                esite = &sites->at(edge->get_end_site_index());

                if(!esite->contains_bwd_edge(edge,true))
                {
                    cout<<"site "<<i<<" has fwd edge from "<<edge->get_start_site_index()<<" to "
                            <<edge->get_end_site_index()<<" but no return\n";
                }
            }
        }

        if( ssite->has_bwd_edge() )
        {
            Edge *edge = ssite->get_first_bwd_edge();
            Site *esite = &sites->at(edge->get_start_site_index());

            if(!esite->contains_fwd_edge(edge,true))
            {
                cout<<"site "<<i<<" has bwd edge from "<<edge->get_start_site_index()<<" to "
                        <<edge->get_end_site_index()<<" but no return\n";
            }

            while( ssite->has_next_bwd_edge() )
            {
                edge = ssite->get_next_bwd_edge();
                esite = &sites->at(edge->get_start_site_index());

                if(!esite->contains_fwd_edge(edge,true))
                {
                    cout<<"site "<<i<<" has bwd edge from "<<edge->get_start_site_index()<<" to "
                            <<edge->get_end_site_index()<<" but no return\n";
                }
            }
        }

    }
}
